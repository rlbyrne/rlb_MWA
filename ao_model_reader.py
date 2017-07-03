#!/usr/bin/python

import sys

def read_model(name):

	filepath = '/nfs/eor-00/h1/rbyrne/MWA/IDL_code/anoko/mwa-reduce/models/model-'+name+'.txt'
	datafile = open(filepath, 'r')
	file_contents = datafile.readlines()
	datafile.close()
	
	source_start_lines = []
	component_start_lines = []
	for i, line in enumerate(file_contents):
		if line.find('source {') != -1:
			source_start_lines.append(i)
		if line.find('component {') != -1:
			component_start_lines.append(i)
	source_start_lines.append(len(file_contents))
	component_start_lines.append(len(file_contents))
	
	source_models = []
	for i in range(len(source_start_lines)-1):
		if file_contents[source_start_lines[i]+1].find('name') != -1:
			name_data = file_contents[source_start_lines[i]+1].split()
			source_name = name_data[1]
		else:
			source_name = ''
		components = []
		for j in range(len(component_start_lines)-1):
			if component_start_lines[j] > source_start_lines[i] and component_start_lines[j] < source_start_lines[i+1]:
				component_data = file_contents[component_start_lines[j]+1: \
							min(component_start_lines[j+1],source_start_lines[i+1])]
				type_name = ''
				ra_deg, dec_deg = [0,0]
				freq_MHz = 0
				flux_I_Jy, flux_Q_Jy, flux_U_Jy, flux_V_Jy = [0,0,0,0]
				for line in component_data:
					if line.find('type') != -1:
						type_data = line.split()
						type_name = type_data[1]
					if line.find('position') != -1:
						pos_data = line.split()
						ra = pos_data[1]
						dec = pos_data[2]
						ra_deg = convert_to_deg(ra)
						dec_deg = convert_to_deg(dec)
					if line.find('frequency') != -1:
						freq_data = line.split()
						freq_MHz = float(freq_data[1])
					if line.find('fluxdensity') != -1:
						flux_data = line.split()
						flux_I_Jy, flux_Q_Jy, flux_U_Jy, flux_V_Jy = [float(flux) for flux in flux_data[2:6]]
				comp = {'type': type_name, 'ra_deg': ra_deg, 'dec_deg': dec_deg, 'freq_MHz': freq_MHz, \
					'flux_I_Jy': flux_I_Jy, 'flux_Q_Jy': flux_Q_Jy, 'flux_U_Jy': flux_U_Jy, 'flux_V_Jy': flux_V_Jy}
				components.append(comp)
				
		source = {'name': source_name, 'components': components}
		source_models.append(source)
		
	return source_models
		
def convert_to_deg(pos):

	#find data type:
	if pos.find('h') != -1:
		hours = True
		hd, minsec = pos.split('h')
	else:
		hours = False
	if pos.find('d') != -1:
		degrees = True
		hd, minsec = pos.split('d')
	else:
		degrees = False
	if hours == degrees:
		print "ERROR: RA/Dec data type is unclear."
		sys.exit(1)
		
	minutes, sec = minsec.split('m')
	hd = float(hd)
	minutes = float(minutes)
	sec = float((sec.split('s'))[0])
	if hd < 0:
		minutes = -minutes
		sec = -sec
	pos_deg = hd + minutes/60. + sec/3600.
	if hours:
		pos_deg = pos_deg/24.*360
	return pos_deg
