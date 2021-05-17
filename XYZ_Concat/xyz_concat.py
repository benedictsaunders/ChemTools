import argparse

#syntax: $ python xyz_concat.py file1 file2 file3 etc...

parser = argparse.ArgumentParser()
parser.add_argument('file', type=argparse.FileType('r'), nargs='+')
args = parser.parse_args()


file_list = []
final_list = []

for arg_f in args.file:
    lines = []
    for line in arg_f:
        lines.append(line.rstrip('\n'))
    del lines[0]
    del lines[0]
    lines.remove('*')
    file_list.append(lines)

for ls in file_list:
    for l in ls:
        final_list.append(l)

new_f = open('new.xyz', 'w+')
new_f.write(str(len(final_list)) + '\n\n')
for l in final_list:
    new_f.write(str(l) + '\n')
new_f.write('*')
new_f.close();
print('Files concatenated.')
