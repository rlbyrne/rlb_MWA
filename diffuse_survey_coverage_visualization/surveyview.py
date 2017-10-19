#!/usr/bin/python

# Tools for characterizing a survey


class Observation:

    def __init__(self, info):
        self.obsid = info[0]
        self.lst = float(info[1])
        self.ra = float(info[2])
        self.dec = float(info[3])
        self.az = float(info[4])
        self.el = float(info[5])


def load_survey(obsfile_name):

    obsfile = open(obsfile_name, "r")
    obsinfo = [line.split("\n") for line in obsfile.readlines()]
    obsfile.close()
    obsinfo = [obs[0] for obs in obsinfo]
    obsinfo = obsinfo[1:]  # remove header
    obsinfo = list(set(obsinfo))
    obsinfo.sort()

    observations = []
    for info in obsinfo:
        info = info.split(", ")
        observations.append(Observation(info))

    return observations


def get_pointings(observations):

    # Note that this code likely doesn't work for surveys other than the
    # diffuse survey.

    # Round the declinations to the nearest 6th to clump them in bands
    decs_round = [int(obs.dec/6.)*6 for obs in observations]
    decs_set = list(set(decs_round))
    decs_set.sort()

    Azimuths = [round(obs.az) for obs in observations]
    Elevations = [round(obs.el) for obs in observations]

    dec_pointings_options = [3, 2, 1, 0, -1, -2, -3]
    azimuth_pointings_options = [2, 1, 0, -1, -2]
    for band_index, dec_band in enumerate(decs_set):
        Azimuths_band = []
        Elevations_band = []
        for obs_index in range(len(observations)):
            if decs_round[obs_index] == dec_band:
                Azimuths_band.append(Azimuths[obs_index])
                Elevations_band.append(Elevations[obs_index])
        AzEls_band = zip(Azimuths_band, Elevations_band)
        AzEls_band_set = list(set(AzEls_band))
        Azimuth_band_set = [value[0] for value in AzEls_band_set]
        Elevations_band_set = [value[1] for value in AzEls_band_set]
        Az_zero = Azimuth_band_set[
            Elevations_band_set.index(max(Elevations_band_set))]
        Azimuth_band_set_sort = list(sorted(Azimuth_band_set))
        Az_wrap = 0
        for Az in Azimuth_band_set_sort:
            if Az > Az_zero + 180:
                Az_wrap += 1
        if Az_wrap > 0:
            Azimuth_band_set_sort = (Azimuth_band_set_sort[
                len(Azimuth_band_set_sort)-Az_wrap:]
                + Azimuth_band_set_sort[:len(Azimuth_band_set_sort)-Az_wrap])

        for obs_index in range(len(observations)):
            if decs_round[obs_index] == dec_band:
                observations[obs_index].pointing = '({}, {})'.format(
                    azimuth_pointings_options[
                        Azimuth_band_set_sort.index(Azimuths[obs_index])],
                    dec_pointings_options[band_index]
                    )

    return observations
