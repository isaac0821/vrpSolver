import vrpSolver
import random
import datetime

def runInstance(instance):
	def readInstanceCVRPTW(fileName):
		f = open(fileName, 'r')
		
		line = f.readline()
		fields = str.split(line)
		instanceName = fields[0]
		
		for i in range(3):
			line = f.readline()
			
		line = f.readline()
		fields = str.split(line)
		vehNum = int(fields[0])
		vehCap = int(fields[1])
		
		for i in range(4):
			line = f.readline()
		
		nodes = {}
		for line in f:
			fields = line.split()
			try:
				nodes[int(fields[0])] = {
					'loc': (float(fields[1]), float(fields[2])),
					'demand': float(fields[3]),
					'tStart': float(fields[4]),
					'tEnd': float(fields[5]),
					'serviceTime': float(fields[6])
				}
			except:
				pass
		
		return instanceName, vehNum, vehCap, nodes

	instanceName, vehNum, vehCap, nodes = readInstanceCVRPTW(instance)

	print("Instance Name", instanceName)
	print("Max number of vehicle", vehNum)
	print("Vehicle Capacity", vehCap)
	print(nodes)

	start = datetime.datetime.now()
	res = vrpSolver.bounds.lb.cgCVRPTW(
		nodes, 
		depotID=0, 
		vehNum=vehNum, 
		vehCap=vehCap,
		cutoffTimeSub=10,
		cutoffTimeMas=180)
	end = datetime.datetime.now()
	print("Done!")
	f = open("Result_CVRPTW.txt", 'a')
	f.write("Instance: %s \tCVRPTW: %s\t EarlyBranching: %s\t runtime: %s\n" % (instanceName, res['ColGen'], res['EarlyBranching'], (end - start).total_seconds()))

if __name__ == "__main__":
	instanceList = [
		"C101.txt",
		"C102.txt",
		"C103.txt",
		"C104.txt",
		"C105.txt",
		"C106.txt",
		"C107.txt",
		"C108.txt",
		"C109.txt",
		"C201.txt",
		"C202.txt",
		"C203.txt",
		"C204.txt",
		"C205.txt",
		"C206.txt",
		"C207.txt",
		"C208.txt",
		"R101.txt",
		"R102.txt",
		"R103.txt",
		"R104.txt",
		"R105.txt",
		"R106.txt",
		"R107.txt",
		"R108.txt",
		"R109.txt",
		"R110.txt",
		"R111.txt",
		"R112.txt",
		"R201.txt",
		"R202.txt",
		"R203.txt",
		"R204.txt",
		"R205.txt",
		"R206.txt",
		"R207.txt",
		"R208.txt",
		"R209.txt",
		"R210.txt",
		"R211.txt",
		"RC101.txt",
		"RC102.txt",
		"RC103.txt",
		"RC104.txt",
		"RC105.txt",
		"RC106.txt",
		"RC107.txt",
		"RC108.txt",
		"RC201.txt",
		"RC202.txt",
		"RC203.txt",
		"RC204.txt",
		"RC205.txt",
		"RC206.txt",
		"RC207.txt",
		"RC208.txt"
	]

	for i in instanceList:
		runInstance(i)