# Reconstructed educational analysis; see README.md and notebook for limitations.
from pathlib import Path
import hashlib
import json
from urllib.request import urlopen
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import curve_fit
from astropy.constants import c
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM

plt.rcParams.update({'figure.dpi': 120, 'font.size': 11})
Path('figures').mkdir(exist_ok=True)
C_KMS = c.to_value(u.km/u.s)

DATA_URL = ('https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/'
            'Pantheon%2B_Data/4_DISTANCES_AND_COVAR/Pantheon%2BSH0ES.dat')
EXPECTED_SHA256 = '1cb0fc379ef066afdc2ffd1857681cc478024570d8a3eba284fb645775198cf8'
data_path = Path('data/Pantheon+SH0ES.dat')
data_path.parent.mkdir(exist_ok=True)
if not data_path.exists():
    with urlopen(DATA_URL, timeout=60) as response:
        raw = response.read()
    if hashlib.sha256(raw).hexdigest() != EXPECTED_SHA256:
        raise ValueError('Upstream data changed: review before updating the checksum.')
    data_path.write_bytes(raw)
if hashlib.sha256(data_path.read_bytes()).hexdigest() != EXPECTED_SHA256:
    raise ValueError('Cached data do not match the reviewed file.')

table = pd.read_csv(data_path, sep=r'\s+', comment='#')
columns = ['zHD', 'MU_SH0ES', 'MU_SH0ES_ERR_DIAG']
values = table[columns].to_numpy(dtype=float)
valid = np.isfinite(values).all(axis=1) & (values[:, 0] > 0) & (values[:, 2] > 0)
z, mu_obs, mu_err = values[valid].T
print(f'Used {len(z)} of {len(table)} light-curve entries.')
print(f'Redshift range: {z.min():.5f} to {z.max():.5f}')

def E(z, Omega_m):
    return np.sqrt(Omega_m*(1 + np.asarray(z))**3 + 1 - Omega_m)

def luminosity_distance(z, H0, Omega_m):
    integral, _ = quad(lambda zp: 1/E(zp, Omega_m), 0, float(z))
    return (1 + z)*C_KMS/H0*integral

def mu_theory(z, H0, Omega_m):
    distances = np.array([luminosity_distance(zi, H0, Omega_m)
                          for zi in np.atleast_1d(z)])
    return 5*np.log10(distances) + 25

popt, pcov = curve_fit(
    mu_theory, z, mu_obs, p0=[70, 0.3], sigma=mu_err,
    absolute_sigma=True, bounds=([1, 0.001], [200, 0.999])
)
H0_fit, Omega_m_fit = popt
H0_err, Omega_m_err = np.sqrt(np.diag(pcov))
print(f'H0 = {H0_fit:.3f} +/- {H0_err:.3f} km/s/Mpc')
print(f'Omega_m = {Omega_m_fit:.4f} +/- {Omega_m_err:.4f}')
print('Errors above are conditional diagonal-only fit errors.')

mu_model = mu_theory(z, *popt)
residuals = mu_obs - mu_model
chi2_diag = np.sum((residuals/mu_err)**2)
dof = len(z) - len(popt)
print(f'Diagonal chi-square / degrees of freedom = {chi2_diag:.2f} / {dof}')
print(f'Diagonal reduced chi-square = {chi2_diag/dof:.3f}')

z_grid = np.geomspace(z.min(), z.max(), 350)
fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
ax.errorbar(z, mu_obs, yerr=mu_err, fmt='.', ms=3, alpha=0.25,
            color='#4477AA', label='Pantheon+SH0ES entries')
ax.plot(z_grid, mu_theory(z_grid, *popt), color='#CC3311', lw=2,
        label='Fitted flat Lambda-CDM (diagonal approximation)')
ax.set(xscale='log', xlabel='Hubble-diagram redshift z',
       ylabel='Distance modulus (mag)', title='Supernova Hubble diagram')
ax.legend(fontsize=9)
ax.grid(alpha=0.2)
fig.savefig('figures/hubble_diagram.png')
plt.show()

fig, ax = plt.subplots(figsize=(8, 4.5), constrained_layout=True)
ax.errorbar(z, residuals, yerr=mu_err, fmt='.', ms=2, alpha=0.16,
            color='#4477AA', label='Individual residuals')
edges = np.geomspace(z.min()*0.99, z.max()*1.01, 11)
for i, (left, right) in enumerate(zip(edges[:-1], edges[1:])):
    mask = (z >= left) & (z < right)
    if mask.sum() == 0:
        continue
    weights = 1/mu_err[mask]**2
    ax.errorbar(np.exp(np.average(np.log(z[mask]), weights=weights)),
                np.average(residuals[mask], weights=weights),
                yerr=1/np.sqrt(weights.sum()), fmt='o', color='#CC3311', ms=4,
                label='Weighted bin means (diagonal errors)' if i == 0 else None)
ax.axhline(0, color='black', ls='--', lw=1)
ax.set(xscale='log', xlabel='Hubble-diagram redshift z',
       ylabel='Observed minus model modulus (mag)', title='Residual diagnostics')
ax.legend(fontsize=9)
ax.grid(alpha=0.2)
fig.savefig('figures/residuals.png')
plt.show()

def age_of_universe(H0, Omega_m):
    integral, _ = quad(lambda zp: 1/((1+zp)*E(zp, Omega_m)), 0, np.inf)
    hubble_time = (1/(H0*u.km/u.s/u.Mpc)).to_value(u.Gyr)
    return hubble_time*integral

age_fit = age_of_universe(H0_fit, Omega_m_fit)
print(f'Age at fitted parameters: {age_fit:.3f} Gyr')
age_rows = []
for om in [0.2, 0.3, 0.351, 0.4]:
    age = age_of_universe(72.97, om)
    age_rows.append({'H0': 72.97, 'Omega_m': om, 'age_Gyr': age})
    print(f'H0 = 72.97, Omega_m = {om:.3f}: {age:.3f} Gyr')

def mu_fixed_Om(z, H0):
    return mu_theory(z, H0, Omega_m=0.3)

fixed_results = {}
for label, mask in [('all', np.ones(len(z), dtype=bool)),
                    ('low_z', z < 0.1), ('high_z', z >= 0.1)]:
    fitted, covariance = curve_fit(
        mu_fixed_Om, z[mask], mu_obs[mask], p0=[70],
        sigma=mu_err[mask], absolute_sigma=True, bounds=([1], [200])
    )
    h0 = float(fitted[0])
    error = float(np.sqrt(covariance[0, 0]))
    fixed_results[label] = {'N': int(mask.sum()), 'H0': h0,
                            'H0_error_diagonal': error}
    print(f'{label:7s}: N={mask.sum():4d}, H0={h0:.3f} +/- {error:.3f} km/s/Mpc')
print('All fits in this section fix Omega_m = 0.3.')

test_z = np.array([0.002, 0.03, 0.1, 0.5, 1.0, 2.2])
reference = FlatLambdaCDM(H0=H0_fit, Om0=Omega_m_fit, Tcmb0=0)
np.testing.assert_allclose(mu_theory(test_z, *popt),
                           reference.distmod(test_z).value, atol=1e-7, rtol=0)
for om in [0.2, 0.3, 0.351, 0.4]:
    factor = 2/(3*np.sqrt(1-om))*np.arcsinh(np.sqrt((1-om)/om))
    analytic_age = factor*(1/(72.97*u.km/u.s/u.Mpc)).to_value(u.Gyr)
    np.testing.assert_allclose(age_of_universe(72.97, om), analytic_age, rtol=1e-9)
for label, mask in [('all', np.ones(len(z), dtype=bool)),
                    ('low_z', z < 0.1), ('high_z', z >= 0.1)]:
    weights = 1/mu_err[mask]**2
    offset = np.average(mu_obs[mask] - mu_fixed_Om(z[mask], 70), weights=weights)
    exact_h0 = 70*10**(-offset/5)
    exact_error = exact_h0*np.log(10)/5/np.sqrt(weights.sum())
    np.testing.assert_allclose(fixed_results[label]['H0'], exact_h0, rtol=1e-6)
    np.testing.assert_allclose(fixed_results[label]['H0_error_diagonal'],
                               exact_error, rtol=1e-5)
assert np.all(np.isfinite(pcov))
print('Passed: Astropy distances, analytic ages, and exact fixed-density fits.')

results = {
    'method': 'Educational diagonal-only, single-redshift approximation',
    'data_sha256': EXPECTED_SHA256, 'N': int(len(z)),
    'joint_fit': {'H0': float(H0_fit), 'H0_error_diagonal': float(H0_err),
                  'Omega_m': float(Omega_m_fit),
                  'Omega_m_error_diagonal': float(Omega_m_err),
                  'age_Gyr': float(age_fit)},
    'fixed_Omega_m_0.3': fixed_results, 'ages_at_H0_72.97': age_rows,
    'chi2_diagonal': float(chi2_diag), 'dof': int(dof),
}
Path('results.json').write_text(json.dumps(results, indent=2) + '\n')
