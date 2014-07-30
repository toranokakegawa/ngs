import numpy as np
import matplotlib.pyplot as plt







def plot_histo(cS1,cS2,cS3,cS4,cS5,cS6,cS7,cS8,cS9,cS10):
	fig = plt.figure()
	ax = fig.add_subplot(111)
	
	## the data
	N = 5
	spadeSize = [cS1,cS2,cS3,cS4,cS5]
	#menStd =   [2, 3, 4, 1, 2]
	raySize = [cS6,cS7,cS8,cS9,cS10]
	#womenStd =   [3, 5, 2, 3, 3]
	
	#max of sizes in order to have the right Y axe
	ymax = max(cS1,cS2,cS3,cS4,cS5,cS6,cS7,cS8,cS9,cS10)
	
	## necessary variables
	ind = np.arange(N)                # the x locations for the groups
	width = 0.35                      # the width of the bars
	
	## the bars
	rects1 = ax.bar(ind, spadeSize, width,
	                color='red')

	rects2 = ax.bar(ind+width, raySize, width,
                    color='black')
	# axes and labels
	ax.set_xlim(-width,len(ind)+width)
	ax.set_ylim(0,ymax+50)
	ax.set_ylabel('Nomber of contigs')
	ax.set_title('Contigs distribution regarding their sizes')
	#xTickMarks = ['Group'+str(i) for i in range(1,6)]
	xTickMarks = ['[0-100]','[100-250]','[250-500]','[500-1000]','[1000-...]']
	ax.set_xticks(ind+width)
	xtickNames = ax.set_xticklabels(xTickMarks)
	plt.setp(xtickNames, rotation=45, fontsize=10)
	
	## add a legend
	ax.legend( (rects1[0], rects2[0]), ('Spades', 'Ray') )
	
	plt.show()