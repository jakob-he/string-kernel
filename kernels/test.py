import time
import gappy_kernel as gk
import numpy as np
import matplotlib.pyplot as plt
from sklearn import svm
from Bio.Seq import Seq
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.metrics import average_precision_score
from sklearn.metrics import precision_recall_curve, precision_score

# Reads in files
def read(fname):
    sequences=[]
    with open(fname,'r') as f:
        for line in f:
            if line[0]!='>':
                sequences.append(line.split()[0])
    return sequences

# Test how fast the matrix is constructed
def speedMatrix(pos,neg,k,g,got_time=False):
    start = time.time()
    gk.extract_spectrum(pos,k=k,g=g)
    gk.extract_spectrum(neg,k=k,g=g)
    print ("Calculated {}-gappypair in {} seconds".format(k, time.time() - start))
    if got_time:
        start = time.time()
        gk.extract_spectrum(pos,k=k,g=g,reverse=True)
        gk.extract_spectrum(neg,k=k,g=g,reverse=True)
        print ("Calculated {}-gappypair with reverse in {} seconds".format(k, time.time() - start))
        start = time.time()
        gk.extract_spectrum(pos,k=k,g=g,include_flanking=True)
        gk.extract_spectrum(neg,k=k,g=g,include_flanking=True)
        print ("Calculated {}-gappypair with flanking in {} seconds".format(k, time.time() - start))
        start = time.time()
        gk.extract_spectrum(pos,k=k,g=g,include_flanking=True,reverse=True)
        gk.extract_spectrum(neg,k=k,g=g,include_flanking=True,reverse=True)
        print ("Calculated {}-gappypair with flanking and reverse in {} seconds".format(k, time.time() - start))

# Test kernel
def trainSVM(pos,neg,k,g,fname,fname2):
    '''cut1=int(0.8*len(pos))
    cut2=int(0.8*len(neg))
    trainX=pos[:cut1]+neg[:cut2]
    trainY=np.append(np.ones(cut1),np.ones(cut2)*(-1))
    testX=pos[cut1:] + neg[cut2:]
    testY=np.append(np.ones(len(pos[cut1:])),np.ones(len(neg[cut2:]))*(-1))
    clf = svm.SVC(C=0.1, kernel='linear', probability=True)
    X=gk.extract_spectrum(trainX,k=k,g=g)
    clf.fit(X,trainY)'''
    spectrum_pos = gk.extract_spectrum(pos, k, g, sparse=False, include_flanking=False)
    spectrum_neg = gk.extract_spectrum(neg, k, g, sparse=False, include_flanking=False)
    print(spectrum_neg,spectrum_pos)
    y = np.concatenate((np.ones(spectrum_pos.shape[0]), -np.ones(spectrum_neg.shape[0])))
    X = np.concatenate((spectrum_pos, spectrum_neg), axis=0)
    X_train, X_test, y_train, y_test = train_test_split(X,y,test_size=0.1,random_state=42,stratify=y)
    print (X_train.shape, X_test.shape)
    print (y_train.shape, y_test.shape)
    start = time.time()
    clf = SVC(C=0.1, kernel='linear', probability=True)
    clf.fit(X_train, y_train)
    print ("Trained linear SVM on {}-spectrum in {} seconds".format(k, time.time() - start))
    y_score = clf.predict_proba(X_test)
    roc_auc = roc_auc_score(y_score=y_score[:,1], y_true=y_test)
    tpr, fpr, _ = roc_curve(y_score=y_score[:,1], y_true=y_test)
    fig = plt.figure(figsize=(14, 8))
    plt.plot(tpr, fpr, label='Performance Spectrum Kernel (AUC: {:0.2f})'.format(roc_auc),
             lw=4, color='orange')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.plot([0, 1], [0, 1], color='navy', lw=4, linestyle='--', label='random')
    plt.legend(loc='lower right')
    fig.savefig(fname)
    print ("Test Set ROC-AUC: {}".format(roc_auc))
    fig = plt.figure(figsize=(14, 8))
    pr_auc = average_precision_score(y_score=y_score[:,1], y_true=y_test)
    pr, rc, _ = precision_recall_curve(y_test, y_score[:,1])
    fig = plt.figure(figsize=(14, 8))
    plt.plot(rc, pr, label='Performance Spectrum Kernel (AUC: {:0.2f})'.format(pr_auc),lw=4, color='orange')
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.plot([0, 1], [0.5, 0.5], color='navy', lw=4, linestyle='--', label='Random')
    plt.ylim([0, 1.05])
    plt.legend(loc='upper right')
    print ("Test Set PR-AUC: {}".format(pr_auc))
    fig.savefig(fname2)

def main():
    posIGF=[Seq(x) for x in read('./testdata/positive_IGF2BP123.fasta')]
    negIGF=[Seq(x) for x in read('./testdata/negative_IGF2BP123.fasta')]
    posPUM=[Seq(x) for x in read('./testdata/positive_PUM2.fasta')]
    negPUM=[Seq(x) for x in read('./testdata/negative_PUM2.fasta')]
    k=3
    g=2
    speedMatrix(posIGF,negIGF,k,g,got_time=True)
    speedMatrix(posPUM,negPUM,k,g)
    trainSVM(posIGF,negIGF,k,g,'./testdata/ROC_AUC_test.pdf','./testdata/PR_AUC_test.pdf')


if __name__ == '__main__':
    main()
