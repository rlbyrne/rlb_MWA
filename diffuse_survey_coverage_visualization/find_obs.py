#!/usr/bin/python

# Script that finds the obsids closest to a given RA/Dec
# Prints both the closest three obsids and the closest three obsids with "best"
# pointings (highest elevation out of observations in its survey with the same
# declination)


import sys
import surveyview


def main(obsfile_name,
         ra_target_h, ra_target_m, ra_target_s,
         dec_target_deg, dec_target_m, dec_target_s):

    ra_target = (int(ra_target_h) + int(ra_target_m)/60.
                 + ra_target_s/3600.)/24.*360.
    dec_target = (int(dec_target_deg) + int(dec_target_m)/60.
                  + dec_target_s/3600.)

    observations = surveyview.load_survey(obsfile_name)
    observations = surveyview.get_pointings(observations)

    for i, obs in enumerate(observations):
        distance_2 = min([
            (obs.dec-dec_target)**2+(use_ra-ra_target)**2 for use_ra in [obs.ra-360, obs.ra, obs.ra+360]
            ])
        obs.distance = distance_2**.5

    obs_sorted = sorted(observations, key=lambda obs: obs.distance)

    print '************'
    print ('CLOSEST OBSERVATIONS TO (RA {:.0f}h{:.0f}m{:.2f}s, '
           'Dec {:.0f}d{:.0f}m{:.2f}s):').format(
                ra_target_h, ra_target_m, ra_target_s,
                dec_target_deg, dec_target_m, dec_target_s)
    print '{:<10}{:>20}'.format('obsid', 'distance (deg)')
    for i in range(3):
        print '{:<10}{:>15.2f}'.format(obs_sorted[i].obsid,
                                       obs_sorted[i].distance)

    for obs in obs_sorted:
        obs.pointing_rank = abs(int(obs.pointing[1:-1].split(', ')[0]))

    print '************'
    print ('CLOSEST OBSERVATIONS TO (RA {:.0f}h{:.0f}m{:.2f}s, '
           'Dec {:.0f}d{:.0f}m{:.2f}s) WITH HIGHEST ELEVATION:').format(
                ra_target_h, ra_target_m, ra_target_s,
                dec_target_deg, dec_target_m, dec_target_s)
    print '{:<10}{:>20}'.format('obsid', 'distance (deg)')
    obs_use = [obs for obs in obs_sorted if obs.pointing_rank == 0]
    for i in range(3):
        print '{:<10}{:>15.2f}'.format(obs_use[i].obsid,
                                       obs_use[i].distance)
    print '************'


if __name__ == '__main__':
    main('/Users/ruby/EoR/sidelobe_survey_obsinfo.txt', 0, 0, 0, -10, 50, 0)
