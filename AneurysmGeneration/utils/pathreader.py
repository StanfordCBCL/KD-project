import xml.etree.ElementTree as ET
import re 

def read_centerline(path_name):

	with open(path_name) as f:
		xml = f.read()
	root = ET.fromstring(re.sub(r"(<\?xml[^>]+\?>)", r"\1<root>", xml) + "</root>")

	for child in root:
		for child_child in child:
			print child_child.tag, child_child.attrib


	path_points = root[1][0][0][1]

	point_list = []
	for point in path_points:
		point_coords = point[0].attrib
		point_list.append( [point_coords['x'], point_coords['y'], point_coords['z']])


	return point_list

if __name__ == "__main__":
	print "testing pathreader.py"
	print "---------------------"

	file_loc = "/Users/alex/Documents/lab/KD-project/AneurysmGeneration/models/SKD0050/"
	file_name = "lad.pth"

	read_centerline(file_loc + file_name)


