import pandas as pd
import matplotlib.pyplot as plt

outDir = "/home/barwu/repos/MuonEffNN/8thTry"

trainDF = pd.read_csv(outDir+"/"+"trainLoss.log", header = None)
validationDF = pd.read_csv(outDir+"/"+"validationLoss.log", header = None)

trainDF['group']=trainDF.index//100
train_avg=trainDF.groupby(['group']).mean()

ax = train_avg.plot(x = 0, y = 2, label = "Training", color='m')
plt.xlabel("Epoch")
plt.ylabel("Loss")
plt.savefig(outDir+"/"+"trainingCurves_training.png")
#plt.show()
validationDF.plot(x = 0, y = 2, label = "Validation", color='g')
plt.xlabel("Epoch")
plt.ylabel("Loss")
plt.savefig(outDir+"/"+"trainingCurves_validating.png")
plt.show()