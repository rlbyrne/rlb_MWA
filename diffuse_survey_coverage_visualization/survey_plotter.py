#!/usr/bin/python

# Code for plotting the coverage of a survey
# Creates the following plots:
#   plot_azels: Creates a scatter plot of Az/El values represented in the
#       survey
#   plot_radecs_colorcode_decs: Creates a scatter plot of RA/Dec values present
#       in the survey; each declination band has a different color
#   generate_radec_animation: Creates a series of plots showing the RA/Dec
#       values of each observation. These plots can be combined to create an
#       animation of the survey in the order the data was taken.
#   plot_radec_pointings_coverage: Creates a heatmap showing the
#       number of different pointings covering points in the RA/Dec plane
#   plot_radec_obs_coverage: Creates a heatmap showing the
#       number of different observations covering points in the RA/Dec plane
#   radec_reference_for_images: Adds a reference plot to the top of FHD output
#       images that shows the observation's field in relation to the survey
#       field, EoR fields, and the brightest sources (to run on all images in
#       a directory, use radec_reference_for_images_wrapper)


import os
import surveyview
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib import gridspec


# Set color scheme for all plots
# Copy array 3 times to make sure there are enough elements
use_colors = ['black', 'red', 'green', 'magenta', 'cyan', 'yellow', 'blue']*3


def plot_azels(obs_info_file, save_loc):

    # Plots the points in the Az/El plane that are covered by the survey
    # Az/El values are rounded to prevent redundancy
    # Data points are color-coded based on their declination

    save_loc = format_save_loc(save_loc, 'AzEls_plot')
    observations = surveyview.load_survey(obs_info_file)

    # Round the declinations to the nearest 6th to clump them in 7 Dec bands
    decs_round = [int(obs.dec/6.)*6 for obs in observations]
    decs_set = list(set(decs_round))
    decs_set.sort()

    plt.figure()
    for i in range(len(decs_set)):
        use_az = []
        use_el = []
        for index in range(len(observations)):
            if decs_round[index] == decs_set[i]:
                use_az.append(observations[index].az)
                use_el.append(observations[index].el)
        use_azels = zip(use_az, use_el)
        use_azels_set = list(set(use_azels))
        use_az_set = [term[0] for term in use_azels_set]
        use_el_set = [term[1] for term in use_azels_set]
        plt.plot(use_az_set, use_el_set, 'o', markersize=10, mfc=use_colors[i],
                 alpha=1)
    plt.xlabel('Azimuth')
    plt.ylabel('Elevation')
    plt.axis([-10, 350, 60, 91])
    plt.grid(True)
    print 'Saving plot to {}'.format(save_loc)
    plt.savefig(save_loc)
    plt.close()


def plot_radecs_colorcode_decs(obs_info_file, save_loc):

    save_loc = format_save_loc(save_loc, 'radec_coverage_colorcode_plot')
    observations = surveyview.load_survey(obs_info_file)

    for obs in observations:
        if obs.ra > 250:
            obs.ra -= 360

    # Round the declinations to the nearest 6th to clump them in 7 Dec bands
    decs_round = [int(obs.dec/6.)*6 for obs in observations]
    decs_set = list(set(decs_round))
    decs_set.sort()

    plt.figure(figsize=(10, 5))
    for i in range(len(decs_set)):
        use_ras = []
        use_decs = []
        for index in range(len(observations)):
            if decs_round[index] == decs_set[i]:
                use_ras.append(observations[index].ra)
                use_decs.append(observations[index].dec)
        plt.plot(use_ras, use_decs, 'o', markersize=10, mfc=use_colors[i],
                 alpha=0.5)
    plt.xlabel('RA')
    plt.ylabel('Dec')
    plt.axis([-100, 200, -65, 15])
    plt.grid(True)
    print 'Saving plot to {}'.format(save_loc)
    plt.savefig(save_loc)
    plt.close()


def generate_radec_animation(obs_info_file, save_dir):

    if save_dir.endswith('.png'):  # User doesn't get to set file name
        save_dir_split = save_dir.split('/')
        save_dir = '/'.join(save_dir_split[:-1])
    if save_dir.endswith('/'):
        save_dir = save_dir[:-1]

    observations = surveyview.load_survey(obs_info_file)

    for obs in observations:
        if obs.ra > 250:
            obs.ra -= 360

    ras = [obs.ra for obs in observations]
    decs = [obs.dec for obs in observations]

    for i, obs in enumerate(obsids):

        if i+1 < 10:
            filepath_num = '00' + str(i+1)
        else:
            if i+1 < 100:
                filepath_num = '0' + str(i+1)
            else:
                filepath_num = str(i+1)
        save_loc = '{}/radec_plot{}.png'.format(save_dir, filepath_num)

        plt.figure(figsize=(17, 5))
        plt.plot(ras[0:i+1], decs[0:i+1], 'o', markersize=90, mfc='blue',
                 alpha=0.03)
        plt.plot(ras[i], decs[i], 'o', markersize=90, mfc='none', mec='red')
        plt.plot(ras[i], decs[i], 'x', markersize=10, mfc='red')
        plt.xticks(range(-100, 200, 10))
        plt.xlabel('RA')
        plt.ylabel('Dec')
        plt.axis('equal')
        plt.axis([-100, 200, -65, 15])
        plt.grid(which='both')
        plt.text(-80, 17, 'Obsid: ' + str(obs))
        print 'Saving plot to {}'.format(save_loc)
        plt.savefig(save_loc)
        plt.close()


def plot_radec_pointings_coverage(obs_info_file, save_loc):

    def check_obsids(observations, ra_target, dec_target, radius,
                     colorbar_max):
        number_pointings = 0
        for i, obs in enumerate(observations):
            distance_2 = (obs.dec-dec_target)**2+(obs.ra-ra_target)**2
            if distance_2 <= radius**2:
                pointings.append(obs.pointing)
        number_pointings = len(set(pointings))
        if number_pointings > colorbar_max:
            number_pointings = colorbar_max
        return number_pointings

    save_loc = format_save_loc(save_loc, 'pointings_coverage')

    ra_range = [-60, 160]
    dec_range = [-70, 20]
    resolution = 1
    radius = 12
    colorbar_max = 10

    ra_vals = [ra_range[0]+resolution*i for i in
               range(int((ra_range[1]-ra_range[0])/resolution))]
    dec_vals = [dec_range[0]+resolution*i for i in
                range(int((dec_range[1]-dec_range[0])/resolution))]
    counts = [[0 for i in range(len(ra_vals))] for i in range(len(dec_vals))]

    observations = surveyview.load_survey(obs_info_file)
    observations = surveyview.get_pointings(observations)

    for i, ra_target in enumerate(ra_vals):
        for j, dec_target in enumerate(dec_vals):
            (counts[-(j+1)])[-(i+1)] = check_obsids(
                observations, ra_target, dec_target, radius, colorbar_max)

    plt.figure(figsize=(9, 3))
    plt.imshow(
        counts,
        interpolation='none',
        extent=[ra_range[1]/360.*24., ra_range[0]/360.*24., dec_range[0],
                dec_range[1]],
        aspect=24/360.,
        cmap=plt.get_cmap('Blues')
        )
    plt.xlabel('Right Ascension (hours)')
    plt.ylabel('Declination (degrees)')
    plt.tight_layout()
    cbar = plt.colorbar(extend='max')
    cbar.set_label('Number of Unique Pointings')
    print 'Saving plot to {}'.format(save_loc)
    plt.savefig(save_loc)
    plt.close()


def plot_radec_obs_coverage(obs_info_file, save_loc):

    def check_obsids(observations, ra_target, dec_target, radius,
                     colorbar_max):
        number_pointings = 0
        for i, obs in enumerate(observations):
            distance_2 = (obs.dec-dec_target)**2+(obs.ra-ra_target)**2
            if distance_2 <= radius**2:
                number_pointings += 1
        if number_pointings > colorbar_max:
            number_pointings = colorbar_max
        return number_pointings

    save_loc = format_save_loc(save_loc, 'obs_coverage')

    ra_range = [-60, 160]
    dec_range = [-70, 20]
    resolution = 1
    radius = 12
    colorbar_max = 15

    ra_vals = [ra_range[0]+resolution*i for i in
               range(int((ra_range[1]-ra_range[0])/resolution))]
    dec_vals = [dec_range[0]+resolution*i for i in
                range(int((dec_range[1]-dec_range[0])/resolution))]
    counts = [[0 for i in range(len(ra_vals))] for i in range(len(dec_vals))]

    observations = surveyview.load_survey(obs_info_file)
    observations = surveyview.get_pointings(observations)

    for i, ra_target in enumerate(ra_vals):
        for j, dec_target in enumerate(dec_vals):
            (counts[-(j+1)])[-(i+1)] = check_obsids(
                observations, ra_target, dec_target, radius, colorbar_max)

    plt.figure(figsize=(9, 3))
    plt.imshow(
        counts,
        interpolation='none',
        extent=[ra_range[1]/360.*24., ra_range[0]/360.*24., dec_range[0],
                dec_range[1]],
        aspect=24/360.,
        cmap=plt.get_cmap('Purples')
        )
    plt.xlabel('Right Ascension (hours)')
    plt.ylabel('Declination (degrees)')
    plt.tight_layout()
    cbar = plt.colorbar(extend='max')
    cbar.set_label('Number of Unique Observations')
    print 'Saving plot to {}'.format(save_loc)
    plt.savefig(save_loc)
    plt.close()


def radec_reference_for_images(obs_info_file, image_filename, save_loc,
                               reduced_obslist=''):

    save_loc = format_save_loc(save_loc, 'radec_reference')

    observations = surveyview.load_survey(obs_info_file)
    for obs in observations:
        if obs.ra > 250:
            obs.ra -= 360

    if reduced_obslist != '':
        reduced_obslist = open(reduced_obslist, "r")
        obs_reduced = reduced_obslist.readlines()
        reduced_obslist.close()
        obs_reduced = [obs.strip() for obs in obs_reduced]
        observations = [
            obs for obs in observations if obs.obsid in obs_reduced]

    ras = [obs.ra for obs in observations]
    decs = [obs.dec for obs in observations]

    highlight_obs = (image_filename.split('/')[-1]).split('_')[0]
    i = [obs.obsid for obs in observations].index(highlight_obs)

    a_teams = surveyview.get_a_team_sources()
    for source in a_teams:
        if source.ra > 250:
            source.ra -= 360

    plt.figure(figsize=(11, 14))
    gs = gridspec.GridSpec(2, 3, height_ratios=[1, 4], width_ratios=[1, 50, 1])

    plt.subplot(gs[0, 1])
    plt.plot(ras, decs, 'o', markersize=90, mfc='0.9', mec='none', zorder=1)
    plt.plot(ras[i], decs[i], 'o', markersize=90, mfc='none', mec='red',
             zorder=11)
    plt.plot(ras[i], decs[i], 'x', markersize=10, mec='red', zorder=12)
    plt.plot([source.ra for source in a_teams],
             [source.dec for source in a_teams],
             'o', markersize=5, mfc='orange', mec='none')
    for index in range(len(a_teams)):
        plt.annotate(
            a_teams[index].name, (a_teams[index].ra, a_teams[index].dec))
    plt.plot([0, 60], [-27, -27], 'o', markersize=90, mfc='cyan', alpha=0.2,
             mec='none', zorder=2)
    plt.annotate('EoR-0', (0, -27))
    plt.annotate('EoR-1', (60, -27))
    plt.xticks(range(-100, 200, 10))
    plt.xlabel('RA')
    plt.ylabel('Dec')
    plt.axis('equal')
    plt.axis([200, -100, -65, 15])
    plt.grid(which='both', zorder=10)

    plt.subplot(gs[1, :])
    img = mpimg.imread(image_filename)
    plt.imshow(img)

    plt.tight_layout()
    print 'Saving plot to {}'.format(save_loc)
    plt.savefig(save_loc, dpi=300)
    plt.close()


def radec_reference_for_images_wrapper(obs_info_file, images_dir, save_dir,
                                       reduced_obslist=''):

    if save_dir.endswith('.png'):  # User doesn't get to set file name
        save_dir_split = save_dir.split('/')
        save_dir = '/'.join(save_dir_split[:-1])
    if save_dir.endswith('/'):
        save_dir = save_dir[:-1]

    if images_dir.endswith('/'):
        images_dir = images_dir[:-1]

    images = os.listdir(images_dir)
    for filename in images:
        if not filename.endswith('.png'):
            images.remove(filename)

    for filename in images:
        output_path = '{}/{}_radec_ref.png'.format(save_dir, filename[0:-4])
        radec_reference_for_images(
            obs_info_file, '{}/{}'.format(images_dir, filename),
            output_path, reduced_obslist)


def format_save_loc(save_loc, default_filename):
    if not save_loc.endswith('.png'):
        if save_loc.endswith('/'):
            save_loc = save_loc[:-1]
        save_loc = '{}/{}.png'.format(save_loc, default_filename)
    return save_loc


if __name__ == '__main__':
    radec_reference_for_images_wrapper('/Users/ruby/EoR/sidelobe_survey_obsinfo.txt', '/Users/ruby/EoR/aws_plots/restored_I', '/Users/ruby/EoR/radec_reference_plots', '/Users/ruby/EoR/diffuse_survey_good_pointings.txt')
