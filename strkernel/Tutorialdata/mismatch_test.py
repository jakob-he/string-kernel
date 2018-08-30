"""
Test for mismatch string kernel
Created on July 10, 2018
@author: Meng Zhang
"""

from mismatch_kernel import MismatchKernel
from mismatch_kernel import preprocess

from Bio import SeqIO
from Bio.Seq import Seq
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score
from sklearn.metrics import classification_report  # classfication summary
import matplotlib.pyplot as plt
import numpy as np
from numpy import random
from timeit import default_timer as timer


# load the data
posSeq = [seq.seq for seq in SeqIO.parse('/notebook_data/positive_IGF2BP123.fasta', 'fasta')]
negSeq = [seq.seq for seq in SeqIO.parse('/notebook_data/negative_IGF2BP123.fasta', 'fasta')]

# a quick check by randomly opting 1000 positive sequences 
# and 1000 negative sequneces to save running time
posX = preprocess(random.choice(posSeq, 1000))
negX = preprocess(random.choice(negSeq, 1000))

# label positive sequences as 1 and negative as 0
posY = np.ones(len(posX), dtype=int)
negY = np.zeros(len(negX), dtype=int)

start = timer()

# (5,1)-mismatch kernel
posKernels = MismatchKernel(l=4, k=5, m=1).get_kernel(posX).kernel
negKernels = MismatchKernel(l=4, k=5, m=1).get_kernel(negX).kernel

end = timer()
print("Time used to compute kernels:", end-start)

# merge data
X = np.concatenate([posKernels,negKernels])
y = np.concatenate([posY,negY])

# split training and test data, choose 30% as test data. 
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.30)

clf = SVC()
clf.fit(X_train, y_train)

y_true, y_pred = y_test, clf.predict(X_test)
print(classification_report(y_true, y_pred))

y_score = clf.decision_function(X_test)

'''Plots a roc curve including a baseline'''
# compute true positive rate and false positive rate
fpr, tpr, thresholds = roc_curve(y_test, y_score)
roc_auc = auc(fpr, tpr)  # compute auc

plt.figure()
lw = 2
plt.plot(fpr, tpr, color='darkorange',
         lw=lw, label='ROC curve (area = %0.2f)' % roc_auc)
plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver operating curve')
plt.legend(loc="lower right")
plt.show()

''' Plot precision and recall curve '''
precision, recall, _ = precision_recall_curve(y_test, y_score)
average_precision = average_precision_score(y_test, y_score)

plt.figure()
plt.step(recall, precision, color='b', alpha=0.2, where='post')
plt.fill_between(recall, precision, step='post', alpha=0.2, color='b')
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.ylim([0.0, 1.05])
plt.xlim([0.0, 1.0])
plt.title('2-class Precision-Recall curve: AP={0:0.2f}'.format(average_precision))
plt.show()
