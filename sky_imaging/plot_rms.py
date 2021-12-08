import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('/Users/rubybyrne/rlb_MWA/sky_imaging')
import calculate_empirical_rm

# Sets math mode text to sans serif
params = {'text.usetex': False, 'mathtext.fontset': 'stixsans'}
plt.rcParams.update(params)

eor0_obslist = [
    '1131454296',
    '1131713632',
    '1131710032',
    '1131456096',
    '1131540824',
    '1131539024',
    '1131538904',
    '1131540704',
    '1131455736',
    '1131710392',
    '1131708952',
    '1131457176',
    '1131716512',
    '1131458976',
    '1131712192',
    '1131453936',
    '1131457536',
    '1131537704',
    '1131543584'
]
obslist_full = [
    '1131551744',
    '1130783824',
    '1131562544',
    '1131709912',
    '1130776864',
    '1131461496',
    '1130782264',
    '1131715432',
    '1131733552',
    '1131542624',
    '1130773144',
    '1131461376',
    '1131557144',
    '1131454296',
    '1131731752',
    '1130778664',
    '1131470496',
    '1131559064',
    '1131717232',
    '1131463536',
    '1130773264',
    '1131463416',
    '1131717352',
    '1131713632',
    '1131478056',
    '1131468936',
    '1131468696',
    '1131535424',
    '1131463296',
    '1131465216',
    '1131710032',
    '1130776624',
    '1131456096',
    '1131540824',
    '1131711952',
    '1131459576',
    '1131477936',
    '1131733672',
    '1131564464',
    '1130787784',
    '1131461616',
    '1131558944',
    '1131470616',
    '1131549944',
    '1131553544',
    '1131459696',
    '1130780464',
    '1131726352',
    '1131470736',
    '1131548024',
    '1131710152',
    '1130785864',
    '1131544424',
    '1131542504',
    #'1131733432',
    '1131735232',
    '1131553664',
    '1131724432',
    '1131542744',
    '1131455976',
    '1131719152',
    '1131454416',
    '1130787544',
    '1130776744',
    '1130780224',
    '1131551624',
    '1131722632',
    '1131547904',
    '1131562664',
    '1131550064',
    '1131537104',
    '1131555224',
    '1131467136',
    '1131539024',
    '1131555344',
    '1131546104',
    '1131548144',
    '1131472416',
    '1131558824',
    '1131544304',
    '1130789584',
    '1131476136',
    '1130789344',
    '1131722872',
    '1130785744',
    '1131730072',
    '1131459816',
    '1131564584',
    '1131457776',
    '1131724552',
    '1130787664',
    '1130778424',
    '1131728152',
    '1131722752',
    '1131538904',
    '1131544544',
    '1130778544',
    '1131467016',
    '1131546344',
    '1130789464',
    '1131713512',
    '1131546224',
    '1131474336',
    '1130782144',
    '1131735472',
    '1130775064',
    '1130774824',
    '1131720832',
    '1130774944',
    '1131557264',
    '1130783944',
    '1131472296',
    '1131465096',
    '1131457896',
    '1131555464',
    '1131562424',
    '1131551864',
    '1131540704',
    '1130780344',
    '1131731632',
    '1131468816',
    '1131472536',
    '1130773024',
    '1131474096',
    '1131465336',
    '1131715552',
    '1131458016',
    '1131540944',
    '1131557024',
    '1131731872',
    '1131553424',
    '1131560864',
    '1130784064',
    '1131466896',
    '1130782024',
    '1131560624',
    '1131474216',
    '1131564344',
    '1131729952',
    '1131560744',
    '1130785624',
    '1131709432',
    '1131536624',
    '1131536384',
    '1131711112',
    '1131709192',
    '1131710992',
    '1131453456',
    '1131565304',
    '1131478776',
    '1131566504',
    '1131565184',
    '1131566624',
    '1131566744',
    '1131565064',
    '1131567944',
    '1131478656',
    '1131568544',
    '1131739432',
    '1130788504',
    '1130788264',
    '1131740752',
    '1131455736',
    '1131710392',
    '1131708952',
    '1131457176',
    '1131716512',
    '1131458976',
    '1131712192',
    '1131453936',
    '1131457536',
    '1131537704',
    '1131543584'
]
rm_file = '/Users/rubybyrne/diffuse_survey_rm_empirical_Jul2021.csv'
rm_data = np.genfromtxt(
    rm_file, delimiter=',', dtype=None, names=True
)
rms = np.array([
    rm_data['RM'][np.where(rm_data['ObsID'] == int(obsid))][0] for obsid in obslist_full
])

plt.rcParams.update({'font.size': 12})
plt.figure(figsize=(6,4))
plt.hist(rms, range=[np.min(rms),0], bins=20, linewidth=0)
plt.xlabel('RM (rad/m$^{2}$)')
plt.ylabel('Histogram count')
plt.tight_layout()
plt.savefig('/Users/rubybyrne/Downloads/RM_hist', dpi=600)

rm_orig_file = '/Users/rubybyrne/diffuse_survey_rm_tot.csv'
rm_orig_data = np.genfromtxt(
    rm_orig_file, delimiter=',', dtype=None, names=True
)
rms_orig = np.array([
    rm_orig_data['RM'][np.where(rm_orig_data['ObsID'] == int(obsid))][0] for obsid in obslist_full
])

rot_angles = calculate_empirical_rm.get_effective_rotation_angles(rms, 167., 198.)
rot_angles_orig = calculate_empirical_rm.get_effective_rotation_angles(rms_orig, 167., 198.)

plt.rcParams.update({'font.size': 12})
plt.figure(figsize=(6,4))
plt.hist(rot_angles, range=[-np.pi, np.pi], bins=20, linewidth=0)
plt.xlim([-np.pi, np.pi])
plt.xlabel('Rotation angle (rad)')
plt.ylabel('Histogram count')
plt.tight_layout()
plt.savefig('/Users/rubybyrne/Downloads/rot_angle_hist', dpi=600)

rm_diff = rms-rms_orig
rm_diff = rm_diff[np.where(rm_diff != 0)]
plt.rcParams.update({'font.size': 12})
plt.figure(figsize=(6,4))
plt.hist(rm_diff, range=[-.3, .3], bins=20, linewidth=0)
plt.xlabel('RM correction (rad/m$^{2}$)')
plt.ylabel('Histogram count')
plt.tight_layout()
plt.savefig('/Users/rubybyrne/Downloads/RM_hist_diff', dpi=600)

rot_angles_diff = rot_angles-rot_angles_orig
rot_angles_diff = rot_angles_diff[np.where(rot_angles_diff != 0)]
rot_angles_diff[np.where(rot_angles_diff < -np.pi)] += 2*np.pi
rot_angles_diff[np.where(rot_angles_diff > np.pi)] -= 2*np.pi
plt.rcParams.update({'font.size': 12})
plt.figure(figsize=(6,4))
plt.hist(rot_angles_diff, range=[-np.pi, np.pi], bins=20, linewidth=0)
plt.xlim([-np.pi, np.pi])
plt.xlabel('Rotation angle correction (rad)')
plt.ylabel('Histogram count')
plt.tight_layout()
plt.savefig('/Users/rubybyrne/Downloads/rot_angle_hist_diff', dpi=600)
