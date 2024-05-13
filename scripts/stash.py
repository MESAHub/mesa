#!/usr/bin/python
from datetime import datetime
from shutil import copytree, copyfile, ignore_patterns
from glob import glob
from os.path import isfile, getmtime
from os import remove, removedirs
import argparse

parser = argparse.ArgumentParser(description='Store a MESA experiment.')
parser.add_argument('-p', metavar='p', type=str, default=None, help='Specifies the starting photo.')
parser.add_argument('-N', metavar='N', type=int, default=5, help='The first N and last N photos by modified-time will be stashed.')
parser.add_argument('-m', metavar='m', type=str, default='', help='Notes about the experiment.')
parser.add_argument('--png', action='store_true', help='Stash png files. If specified png files will be stashed and then removed from the working directory.')
args = parser.parse_args()

# Name the experiment
now = datetime.now()
time = now.strftime('%Y_%m_%d_%H_%M_%S')
experiment_name = 'experiment_' + time
experiment_folder = 'experiments/'+experiment_name

copied_pngs = []

# Stash everything that isn't a png, movie, photo, or compiled binary file.
def ignore(path, names):
	global copied_pngs

	if './experiments' in path or './photos' in path or '__pycachce__' in path:
		return names
	ig = set()

	if not args.png:
		# Ignore png's
		ig.update(n for n in names if n[-4:] == '.png')
	else:
		# Note which png's were copied so we can remove them later.
		copied_pngs = copied_pngs + list(path + '/' + n for n in names if n[-4:] == '.png')

	ig.update(n for n in names if n[-4:] == '.mov')
	ig.update(n for n in names if n[-2:] == '.o')
	ig.update(n for n in names if n[-4:] == '.mod' and path == './make')
	ig.update(n for n in names if n[-5:] == '.smod')
	ig.update(n for n in names if n == 'star')
	ig.update(n for n in names if n == 'stash.py')
	return ig


copytree('./', experiment_folder, ignore=ignore)

# Stash the first N photos and the last N photos by last-modified time.
N = args.N
photo_dir = glob('photos/*')
photos = list(photo for photo in photo_dir if isfile(photo))
photo_times = list(getmtime(photo) for photo in photos)
photos_and_times = zip(*(photos, photo_times))
sorted_photos = sorted(photos_and_times, key=lambda x: x[1])
sorted_photos = list(photo[0] for photo in sorted_photos)
ending_photo = sorted_photos[-1].lstrip('./photos/')
to_stash = list(set(sorted_photos[:N] + sorted_photos[-N:]))

# Stash the starting photo, if specified
if args.p is not None:
	to_stash.append('photos/' + args.p)
	to_stash = list(set(to_stash))

for photo in to_stash:
	copyfile(photo, experiment_folder + '/' + photo)

# Record experiment notes
with open('experiments/' + experiment_name + '.md','w') as log:
	log.write('## ' + experiment_name + '\n')
	if args.p is not None:
		log.write('- Starting photo: ' + args.p + '\n')
	log.write('- Ending photo: ' + ending_photo + '\n')
	if args.m != '':
		log.write(args.m + '\n')
	log.write('\n')

# Remove contents of LOGS
log_files = glob('LOGS/*')
for f in log_files:
	if isfile(f):
		remove(f)

# If we stashed png's, remove those from the main directory.
for f in copied_pngs:
	remove(f)

# Remove empty experiments directory from the new experiment
removedirs(experiment_folder + '/experiments')

print('Stashed ' + experiment_name + '.')

