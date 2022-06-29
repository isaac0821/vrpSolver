CONST_EPSILON = 0.00001

CONST_EARTH_RADIUS_MILES = 3958.8
CONST_EARTH_RADIUS_METERS = 6378137.0

config = {
	'MESSAGE_SHOW_WARNING': False,
	'MESSAGE_SHOW_ERROR': True,
	'DEFAULT_ASYMFLAG': False
}

def msgWarning(*args):
	if (config['MESSAGE_SHOW_WARNING'] and config['MESSAGE_SHOW_ERROR']):
		print(*args)
	return

def msgError(*args):
	if (config['MESSAGE_SHOW_ERROR']):
		print(*args)
	return
