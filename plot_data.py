# -------------------------------------------------------
#
#   PySPT, a likelihood framework for fitting PMF model to
#       South Pole Telescope with Python.
#
#   file: emc_spt.py
#
#       This file contains the main running code for the SPT python
#       project.
#
#
#   @author: Alex Zucca, azucca@sfu.ca, azucca@dwavesys.com
#
#   ------------------
#   MODEL PARAMETERS LIST:
#
#   0. Abb:        lensing amplitude
#   1. r:          tensor to scalar ratio
#   2. constant:   constant to add to the cls
#   3. Add:            dust amplitude,
#   4. Poisson150_bf:  poisson noise at 150GHz amplitude
#   5. Poisson_90x150: poisson noise at 90x150GHz amplitude
#   6. Poisson90:      poisson noise at 90GHz amplitude,
#   7. MapBcal150:     map calibration amplitude at 150GHz
#   8. MapBcal90:      map calibration amplitude at 90GHz,
#   9. bf1:            beam factor 1
#   :   ...         ...
#   15. bf7:            beam factor 7
#   16. Bpmf:           PMF amplitude
#   17. beta_pmf        ratio for generation time
#
#
# -------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import corner
import emcee
import python_spt.pyspt as spt


##########################
# SET PARAMETERS
helical = False
magind = 2

# set the random seed
np.random.seed(0)

# fix a few parameters
spt_l_min = 3
spt_l_max = 3000

# MCMC parameters
mcmc_steps = 25000

# n_walkers should be at least twice the num of params
n_walkers = 40

# parameters to plot, follow order above
params_to_plot = [0, 1, 2, 16, 17]

# the following is automatically done
labels = [r'$A_{\rm lens}$', r'$r$', r'$C$', r'$A_{\rm dust}$', r'$A_P^{150}$', r'$A_P^{90 \times 150}$', r'$A_P^{90}$',
          r'$K^{150}$', r'$K^{90}$', r'$b_f^1$', r'$b_f^2$', r'$b_f^3$', r'$b_f^4$', r'$b_f^5$', r'$b_f^6$', r'$b_f^7$',
          r'$B_{1\, {\rm Mpc}}$', r'$\beta$']

if magind == 2:
    labels[16] = r'$\log_{10} B_{1\, {\rm Mpc}}$'

labels_to_plot = [labels[p] for p in params_to_plot]

##########################
# LOAD THE DATA

def load_data(filename):
    """
    Function to read arrays in Fortran contiguous ordering

    :param filename:
    :return:
    """
    print('loading:', filename)
    return np.asfortranarray(np.loadtxt(filename))


cov = load_data('./data/spt_500d_BB_cov.txt')
win = load_data('./data/spt_500d_BB_windows.txt')
bp = load_data('./data/sptpol_500d_BB_bandpowers.txt')
berr = load_data('./data/spt_500d_BB_beam_error.txt')

# get the bandpowers
spec = np.asfortranarray(bp[:, [5, 4, 3]])

print('corner elements of cov')
print(cov[0, 0])
print(cov[20, 0])
print(cov[0, 20])
print(cov[20, 20])


print('cholesky decomposition of cov')
cov = spt.matrix_cholesky_wrapper(cov)
print('done')
print('corner elements of cov')
print(cov[0, 0])
print(cov[20, 0])
print(cov[0, 20])
print(cov[20, 20])


##########################
# LOAD THE TEMPLATES

tensor_primary = load_data('./templates/primary_tensCls.dat')
# lensed_BB = load_data('./templates/primarylensing_lensedCls.dat')
lensed_BB = load_data('./templates/planck_2018_lensedCls.dat')

# load PMF templates
if magind == -2.9:
    root = 'nonhelicalm3ng29'

else:
    root = 'nonhelicalm3p2'

# Load the templates
pmf_template_tens_m2 = load_data('./templates/'+root+'m2_tensCls.dat')
pmf_template_vect_m1 = load_data('./templates/'+root+'m1_vecCls.dat')

if helical:
    # load the helical templates
    root = root[3:]
    pmf_template_tens_hm2 = load_data('./templates/' + root + 'm2_tensCls.dat')
    pmf_template_vect_hm1 = load_data('./templates/' + root + 'm1_vecCls.dat')


######################
# PREPARE THE Cls
spt_slice = spt_l_max-spt_l_min+1
lowl = 1
highl = lowl+spt_slice

dls_tens = tensor_primary[lowl:highl, 3]
dls_lens = lensed_BB[lowl:highl, 3]
dls_v1 = pmf_template_vect_m1[lowl:highl, 3]
dls_t2 = pmf_template_tens_m2[lowl:highl, 3]

dls_hv1 = None
dls_ht2 = None

if helical:
    dls_hv1 = pmf_template_vect_hm1[lowl:highl, 3]
    dls_ht2 = pmf_template_tens_hm2[lowl:highl, 3]


# Prepare the poisson noise and dust foregrounds
ells = np.array([float(l) for l in range(3, 3001)])
dls_poiss = ells*(ells+1.0)/(3000. * 3001)
dls_galdust = (80.0/ells)**0.58

plt.plot(ells, 0.1*dls_poiss, label='poisson')
plt.plot(ells, 0.0094*dls_galdust, label='dust')
plt.plot(ells, dls_tens, label='tensor')
plt.plot(ells, dls_lens, label='lensing')
markers = ['o', 'v', 's']
labels_data = ['95x95', '150x95', '150x150']
colors = ['C4', 'C5', 'C6']
for i in range(3):
    plt.scatter(bp[:, 0], bp[:, 3+i], marker=markers[i], label=labels_data[i], color=colors[i])
plt.legend()
plt.ylim(-0.5, 0.5)
plt.xlim(0, 2500)
plt.show()



def check_params_limit(theta):
    # extract the parameters
    Abb, r, constant, Add, Poisson150, Poisson_90x150, Poisson90, MapBcal150, MapBcal90, bf1, bf2, bf3, bf4, bf5, \
    bf6, bf7, Bpmf, beta_pmf = theta

    within_bounds = True

    # check the parameters limit
    if Abb < 0 or Abb > 100:
        return False

    if r < 0 or r > 1:
        return False

    if beta_pmf > 40 or beta_pmf < 10:
        return False

    if magind == 2:
        if Bpmf<-3 or Bpmf>1:
            return False
    else:
        if Bpmf < 0 or Bpmf > 10:
            return False

    if constant < -0.5 or constant > 0.5:
        return False

    if Poisson_90x150 < 0 or Poisson_90x150 > 10:
        return False

    if Poisson150 < 0 or Poisson150 > 10:
        return False

    if Poisson90 < 0 or Poisson90 > 10:
        return False

    beam_factors = np.array([bf1, bf2, bf3, bf4, bf5, bf6, bf7])

    for bf in beam_factors:
        if bf < -10 or bf > 10:
            return False

    return within_bounds


############################
# DEFINE THE SPT LIKELIHOOD
def spt_log_likelihood(theta, dls_tens, dls_lens, dls_v1, dls_t2, dls_hv1, dls_ht2, dls_poiss, dls_galdust,
                       helical, cov, win, spec, beam_err):

    # extract the parameters
    Abb, r, constant, Add, Poisson150, Poisson_90x150, Poisson90, MapBcal150, MapBcal90, bf1, bf2, bf3, bf4, bf5,\
    bf6, bf7, Bpmf, beta_pmf = theta

    if not check_params_limit(theta):
        return -np.infty

    # start from LCDM
    # rescale the r somehow...
    # dls_lcdm = dls_lens*Abb + dls_tens*r + constant
    dls_lcdm = dls_lens + dls_tens * r + constant
    # dls_lcdm *= 0.0

    # Mix the PMF Dls here
    # rescale the PMF amplitude (in nG)
    if magind == 2:
        Apmf = (10**(Bpmf) / 3.) ** 4
    else:
        Apmf = (Bpmf/3.)**4

    # rescale the generation time of PMFs
    rat = (beta_pmf/17)**(1.9)

    dls_pmf = Apmf*dls_v1 + Apmf*rat*dls_t2
    if helical:
        # changed this:
        dls_pmf = Apmf*dls_hv1 + Apmf*rat*dls_ht2

    # dls_pmf = dls_pmf * 0.0

    # we'll pass these to the fortran code
    poiss_levels = np.array([Poisson150, Poisson_90x150, Poisson90])
    beam_factors = np.array([bf1, bf2, bf3, bf4, bf5, bf6, bf7])

    # plt.plot(np.arange(3, 3001), 2*dls_lcdm)
    # for i in range(3):
    #     plt.scatter(bp[:, 0], bp[:, 3+i])
    # plt.ylim(-0.5, 0.5)
    # plt.xlim(0, 2500)
    # plt.title('Lensing only')
    # plt.show()
    #
    # plt.plot(np.arange(3, 3001), dls_lcdm + r * dls_tens)
    # for i in range(3):
    #     plt.scatter(bp[:, 0], bp[:, 3 + i])
    # plt.ylim(-0.5, 0.5)
    # plt.xlim(0, 2500)
    # plt.title('+ tensor BB')
    # plt.show()

    # call the Fortran Likelihood
    spt_like = -spt.spt_lnlike_wrapper(dls_lcdm, dls_pmf, dls_poiss, dls_galdust, Add, poiss_levels,
                                      MapBcal150, MapBcal90, beam_factors, beam_err, spec, win, cov)

    # print('spt_like:', spt_like)

    return spt_like


############################
# TEST THE LIKELIHOOD
print('Testing SPT likelihood')
theta = 1, 0.1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 17

print('with parameters:', theta)

spt_like = spt_log_likelihood(theta, dls_tens, dls_lens, dls_v1, dls_t2, dls_hv1, dls_ht2, dls_poiss, dls_galdust,
                              helical, cov, win, spec, berr)

print('spt_likelihood:', spt_like)


#############################
# SET UP MCMC
# best fit from paper
# Abb, r, constant, Add, Poisson150, Poisson_90x150, Poisson90, MapBcal150, MapBcal90,
# bf1, bf2, bf3, bf4, bf5, bf6, bf7, Bpmf, beta_pmf

Abb_bf = 1.0
r_bf = 0.01
constant_bf = 0.0
Add_bf = 0.0094
Poisson150_bf = 0.1
Poisson_90x150_bf = 0.1
Poisson90_bf = 0.1
MapBcal150_bf = 1.0
MapBcal90_bf = 1.0
bf_bf = 0.0
if magind == 2:
    Bpmf_bf = -3
else:
    Bpmf_bf = 0
beta_pmf_bf = 17
args_like = (dls_tens, dls_lens, dls_v1, dls_t2, dls_hv1, dls_ht2,
             dls_poiss, dls_galdust, helical, cov, win, spec, berr)

# start around best fit
bf_theta = np.array([Abb_bf, r_bf, constant_bf, Add_bf, Poisson150_bf, Poisson_90x150_bf,
                     Poisson90_bf, MapBcal150_bf, MapBcal90_bf, bf_bf, bf_bf, bf_bf, bf_bf, bf_bf, bf_bf,
                     bf_bf, Bpmf_bf, beta_pmf_bf])

############################
# TEST THE LIKELIHOOD
print('Testing SPT likelihood with BF')
print('with parameters:', bf_theta)

spt_like = spt_log_likelihood(bf_theta, dls_tens, dls_lens, dls_v1, dls_t2, dls_hv1, dls_ht2, dls_poiss, dls_galdust,
                              helical, cov, win, spec, berr)

print('spt_likelihood:', spt_like)


n_par = bf_theta.shape[0]


# create initial points
initial_points = bf_theta + 1e-2 * np.random.randn(n_walkers, n_par)

# initialize the sampler
# sampler = emcee.EnsembleSampler(n_walkers, n_par, spt_log_likelihood, args=args_like,
#                                 moves=[(emcee.moves.DEMove(), 0.25),
#                                        (emcee.moves.DESnookerMove(), 0.25),
#                                        (emcee.moves.WalkMove(), 0.5), ])

sampler = emcee.EnsembleSampler(n_walkers, n_par, spt_log_likelihood, args=args_like)

# run MCMC
sampler.run_mcmc(initial_points, mcmc_steps, progress=True, skip_initial_state_check=True)

# get the samples
# flat_samples = sampler.get_chain(discard=1000, thin=10, flat=True)
flat_samples = sampler.get_chain(discard=1000, thin=10, flat=True)[:, params_to_plot]

# tau = sampler.get_autocorr_time()
# print('max autocorrelation time:', np.max(tau))

# fig = corner.corner(flat_samples, labels=labels_to_plot, quantiles=[0.16, 0.5, 0.84], show_titles=True)
fig = corner.corner(flat_samples, labels=labels_to_plot)
plt.savefig('./test.pdf', bbox_inches='tight')
# plt.show()

for i in range(len(params_to_plot)):
    mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
    q = np.diff(mcmc)
    # txt = "\mathrm{{{3}}} = {0:.3f}_{{-{1:.3f}}}^{{{2:.3f}}}"
    # txt = txt.format(mcmc[1], q[0], q[1], labels[i])
    print('param', params_to_plot[i], ' BF:', mcmc[1], q[0], q[1])