#!/usr/bin/python


import surveyview
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


# Set color scheme for all plots
# Copy array 3 times to make sure there are enough elements
use_colors = ['black', 'red', 'green', 'magenta', 'cyan', 'yellow', 'blue']*3


def plot_azels(obsfile_name, save_loc):

    # Plots the points in the Az/El plane that are covered by the survey
    # Az/El values are rounded to prevent redundancy
    # Data points are color-coded based on their declination

    save_loc = format_save_loc(save_loc, 'AzEls_plot')
    observations = surveyview.load_survey(obsfile_name)

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


def plot_radecs_colorcode_decs(obsfile_name, save_loc):

    save_loc = format_save_loc(save_loc, 'radec_coverage_colorcode_plot')
    observations = surveyview.load_survey(obsfile_name)

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


def generate_radec_animation(obsfile_name, save_loc):

    if save_loc.endswith('.png'):  # User doesn't get to set file name
        save_loc_split = save_loc.split('/')
        save_loc = '/'.join(save_loc_split[:-1])
    if save_loc.endswith('/'):
        save_loc = save_loc[:-1]

    observations = surveyview.load_survey(obsfile_name)

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
        save_filepath = '{}/radec_plot{}.png'.format(save_loc, filepath_num)

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
        print 'Saving plot to {}'.format(save_filepath)
        plt.savefig(save_filepath)
        plt.close()


def plot_radec_pointings_coverage(obsfile_name, save_loc):

    def check_obsids(observations, ra_target, dec_target, radius,
                     colorbar_max):
        pointings = []
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

    observations = surveyview.load_survey(obsfile_name)
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
    print 'Saving plot to {}'.format(save_filepath)
    plt.savefig(save_filepath)
    plt.close()


def plot_radec_obs_coverage(obsfile_name, save_loc):

    def check_obsids(observations, ra_target, dec_target, radius,
                     colorbar_max):
        pointings = []
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

    observations = surveyview.load_survey(obsfile_name)
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
    print 'Saving plot to {}'.format(save_filepath)
    plt.savefig(save_filepath)
    plt.close()


def format_save_loc(save_loc, default_filename):
    if not save_loc.endswith('.png'):
        if save_loc.endswith('/'):
            save_loc = save_loc[:-1]
        save_loc = '{}/{}.png'.format(save_loc, default_filename)
    return save_loc
