import scipy
import numpy as np

def gaus(x,a,x0,sigma):
	return a*np.exp(-(x-x0)**2/(2*sigma**2))

def fitGaussian(x,y):
	from scipy.optimize import curve_fit
	n = len(x)                          #the number of data
	mean = sum(x*y)/n                   #note this correction
	sigma = sum(y*(x-mean)**2)/n        #note this correction
	try:
		popt,pcov = curve_fit(gaus,x,y,p0=[max(y),x[np.argmax(y)],0.12],maxfev=1000000)
		return popt,pcov
	except RuntimeError:
		print "Error - curve_fit failed"

def makeStrainMPI(alldata,roi,bins,length,xr):
	from mpi4py import MPI
	from mpi4py.MPI import ANY_SOURCE
	import time
	start = time.time()
	comm = MPI.COMM_WORLD
	rank = comm.Get_rank()
	size = comm.Get_size()

	print rank,size

	def strainRange(data_part,xr):
		strainpic = np.zeros((np.shape(data_part[0,:,:])))
		for i in range(len(data_part[0,:,0])):
			for j in range(len(data_part[0,0,:])):
				popt,pcov = fitGaussian(xr,data_part[:,i,j])
				if popt[1] > -0.05 and popt[1] < 0.05:
					strainpic[i,j] = popt[1]/(10.992-7.86E-4)
		strainpic[0,0] = rank
		return strainpic

	ypix = (roi[1]-roi[0])/bins
	xpix = (roi[3]-roi[2])/bins
	
	strainpic = np.zeros((xpix,ypix))

	theta = np.arange(-0.0025*(length/2),0.0025*(length/2)+.001,.0001)
	
	pop = []

	local_n = ypix/size
	local_data = alldata[:,:,rank*local_n:(rank+1)*local_n]

	strainpic_part = strainRange(local_data,xr)

	recv_buffer = np.zeros((np.shape(alldata[:,:,rank*local_n:(rank+1)*local_n])))
	strainpic = np.zeros((np.shape(alldata[0,:,:])))

	output = np.zeros((1,len(strainpic_part[:,0]),len(strainpic_part[0,:])))
	#recv_buffer = np.zeros((len(strainpic_part[:,0]),len(strainpic_part[0,:])))

	# output[0,0,0] = rank
	# output[0,:,:] = strainpic_part
	# print output[0,0,0]

	

	if rank == 0:
		datarank = strainpic_part[0,0]
		print datarank
		strainpic_part[0,0] = 0
		strainpic[:,rank*local_n:(rank+1)*local_n] = strainpic_part
		for i in range(1, size):
			comm.Recv(recv_buffer, ANY_SOURCE)
			datarank = recv_buffer[0][0,0]
			print datarank
			recv_buffer[0][0,0] = 0
			strainpic[:,datarank*local_n:(datarank+1)*local_n] = recv_buffer[0]
	else:
		# all other process send their result
		comm.Send(strainpic_part)

	# root process prints results
	if comm.rank == 0:
		#end = time.time()
		#print "Time:",end - start
		print "Sum:",np.mean(strainpic)
		return strainpic
