import vtk
import numpy as np
from vtk.util import numpy_support as nps
import os
import codecs
import subprocess

modelID = 'KDR32'
#SimModFlag = 'Rigid'
SimModFlag = 'FSI'
path = '/home/ngrandeg/Documents/Toronto_KDProject/' + modelID
AllResultsFileName = path +'/SimResults/' + SimModFlag + '/all_results.vtp'
MeshDirectory = path + '/Mesh/FineMeshBL/mesh-complete/mesh-surfaces/'
WallsMeshName = path + '/Mesh/FineMeshBL/mesh-complete/walls_combined.vtp'



#---------------------------------------------------------------------------------------------------------

def Get_TAWSS(NoW, NameList, MeshDirectory, AllResultsFile):

	TNoP = 0
	c=0

	print 'Reading walls ...'

	for i in range(NoW):


		line = NameList[i]
		WallName =  line.strip()
		FileName = MeshDirectory + WallName

		print FileName

		#Read vtp mesh files
		datareader = vtk.vtkXMLPolyDataReader()
		datareader.SetFileName(FileName)
		datareader.Update()
		data = datareader.GetOutput()

		NoP = data.GetNumberOfPoints()
		print NoP
		TNoP = TNoP + NoP


	print 'Total Number of Points:' + str(TNoP)


	total_tawss =np.zeros((1,TNoP))
	total_osi =np.zeros((1,TNoP))
	total_point =np.zeros((TNoP,3))
	total_GlobalID = []

	for i in range(NoW):


		line = NameList[i]
		WallName =  line.strip()
		FileName = MeshDirectory + WallName

		print FileName

		#Read vtp mesh files
		datareader = vtk.vtkXMLPolyDataReader()
		datareader.SetFileName(FileName)
		datareader.Update()
		data = datareader.GetOutput()
		#print data
		NoP = data.GetNumberOfPoints()
		print NoP

		#Read results file
		resultsreader = vtk.vtkXMLPolyDataReader()
		resultsreader.SetFileName(AllResultsFile)
		resultsreader.Update()
		results = resultsreader.GetOutput()


		tawss_np = np.zeros((1,TNoP))
		osi_np = np.zeros((1,TNoP))
		point_np = np.zeros((TNoP,3))

		gNid_array = nps.vtk_to_numpy(results.GetPointData().GetArray('GlobalNodeID'))


		for k in range(NoP):

			gNid=data.GetPointData().GetArray('GlobalNodeID').GetTuple(k)

			total_GlobalID.append( np.int(gNid[0]))

		        index = np.where(gNid_array == np.int(gNid[0]))[0]

 			tawss_np[0,k]=np.double(results.GetPointData().GetArray('vTAWSS').GetTuple(index[0]))
			osi_np[0,k]=np.double(results.GetPointData().GetArray('vOSI').GetTuple(index[0]))
			total_point[c,:] = np.double(data.GetPoints().GetPoint(k))

			total_tawss[0,c] = np.double(tawss_np[0,k])
			total_osi[0,c] = np.double(osi_np[0,k])
			print total_tawss[0,c]
			c=c+1




		tawss_vtk = nps.numpy_to_vtk(tawss_np[0])
		tawss_vtk.SetName('TAWSS')

		osi_vtk = nps.numpy_to_vtk(osi_np[0])
		osi_vtk.SetName('OSI')

		data.GetPointData().AddArray(tawss_vtk)
		data.GetPointData().AddArray(osi_vtk)
		data.GetCellData().RemoveArray('GlobalElementID2')
		data.Update()



		new=vtk.vtkXMLPolyDataWriter()
		new.SetInput(data)
		new.SetFileName('Results_'+ WallName )
		new.Write()

	av_tawss = np.mean(total_tawss)

	return (total_tawss, total_GlobalID, av_tawss)



def Write_TAWSS(meshdata, total_tawss, total_GlobalID):

	NoMP = meshdata.GetNumberOfPoints()

	meshTAWSS =np.zeros((1,NoMP))

	for k in range(NoMP):

		gMNid=meshdata.GetPointData().GetArray('GlobalNodeID').GetTuple(k)

		#print int(gMNid[0])
		check = np.int(gMNid[0]) in total_GlobalID
		#print check
		if check == True:

			index = total_GlobalID.index( np.int(gMNid[0]))
			#print index

			meshTAWSS[0,k] = total_tawss[0,index]
			#print meshTAWSS[0,k]
		else:

			meshTAWSS[0,k] = 0

	return meshTAWSS

#------------------------------------------------------------------------------------------------------



dummy = subprocess.Popen('ls ' + MeshDirectory +' wall*',shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)


dummylines = dummy.stdout.readlines()
NoL = len(dummylines)


dummylines_RCA=filter(lambda x: x[: +8] =='wall_rca' or x[: +14] =='wall_blend_rca' ,dummylines)
N_RCA=len(dummylines_RCA)

dummylines_LCA=filter(lambda x: x[: +8] =='wall_lca',dummylines)
N_LCA=len(dummylines_LCA)

dummylines_br=filter(lambda x: x[: +10] =='wall_aorta',dummylines)
N_br=len(dummylines_br)

dummylines_aorta=filter(lambda x: x[: +14] == 'wall_aorta.vtp',dummylines)

N_aorta=len(dummylines_aorta)

#Compute results right branch

total_tawss_rca, total_GlobalID_rca, av_tawss_rca = Get_TAWSS(N_RCA, dummylines_RCA, MeshDirectory, AllResultsFileName)
total_tawss_lca, total_GlobalID_lca, av_tawss_lca = Get_TAWSS(N_LCA, dummylines_LCA, MeshDirectory, AllResultsFileName)
total_tawss_aorta, total_GlobalID_aorta, av_tawss_aorta = Get_TAWSS(N_aorta, dummylines_aorta, MeshDirectory, AllResultsFileName)
print dummylines_aorta

#Read vtp wall combined mesh files
meshreader = vtk.vtkXMLPolyDataReader()
meshreader.SetFileName(WallsMeshName)
meshreader.Update()
meshdata = meshreader.GetOutput()

meshTAWSS_rca = Write_TAWSS(meshdata, total_tawss_rca, total_GlobalID_rca)
meshTAWSS_lca = Write_TAWSS(meshdata, total_tawss_lca, total_GlobalID_lca)

RCA_tawss_vtk = nps.numpy_to_vtk(meshTAWSS_rca[0])
RCA_tawss_vtk.SetName('TAWSS_rca')
meshdata.GetPointData().AddArray(RCA_tawss_vtk)
meshdata.GetCellData().RemoveArray('GlobalElementID2')
meshdata.Update()
LCA_tawss_vtk = nps.numpy_to_vtk(meshTAWSS_lca[0])
LCA_tawss_vtk.SetName('TAWSS_lca')
meshdata.GetPointData().AddArray(LCA_tawss_vtk)
meshdata.Update()


print 'RCA branch tawss (point average) = ' + str(av_tawss_rca)
print 'LCA branch tawss (point average) = ' + str(av_tawss_lca)


print 'writing vtp file ...'
rca=vtk.vtkXMLPolyDataWriter()
rca.SetInput(meshdata)
rca.SetFileName('LocalResults.vtp')
rca.Write()

#surf= vtk.vtkMassProperties()
#surf.SetInput(meshdata)
#surf.Update()
#area = surf.GetSurfaceArea()
#print area
