

def initialize(): 
	global PTM_CONTACT_DICT 
	#to look up a ligand code replace the <> with intention in this URL
	# https://www.rcsb.org/ligand/<ALY>
	PTM_CONTACT_DICT = {'PTR':'Y', 'MSE':'S', 'SEP':'S', 'CAS':'C', 'ALY':'K'}
	#FYQ, AYI, YEN are synthetic pTYR sequences that act as inhibitors, have to decide what to do about that.
	global MOLECULES_TO_REMOVE
	MOLECULES_TO_REMOVE = ['HOH', 'SO4', 'NBS', 'CSO', 'MN', 'MG', 'ZN', 'MYR', 'P16', 'GOL', 'QUE', 'CA', 'ANP', 'EDO', 'DVT', 'CL', 'PO4', 'FMT', 'ACE', 'CAT', '1N1', 'VSH', 'PB']
	synthetic_inhibitors = ['FYQ', 'AYI', 'YEN', 'AYQ']
	MOLECULES_TO_REMOVE = MOLECULES_TO_REMOVE+synthetic_inhibitors