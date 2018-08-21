#!/usr/bin/env python3
'''
Example application of the motif kernel with a SVM.
'''

# preprocessing
from Bio.Seq import Seq
import Bio.SeqIO as sio
import numpy as np
import time

# SVM
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report  # classfication summary
# ROC and precision-recall curve
from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score

# plotting
import matplotlib.pyplot as plt  # plotting

# own libraries
from kernels.motif_kernel import motifKernel


def plot_roc_curve(y_test, y_score):
    '''Plots a roc curve including a baseline'''
    fpr, tpr, thresholds = roc_curve(y_test, y_score)
    roc_auc = auc(fpr, tpr)

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


def plot_prec_recall_curve(y_test, y_scores):
    '''Plots a precision-recall curve including a baseline'''
    precision, recall, thresholds = precision_recall_curve(y_test, y_scores)
    average_precision = average_precision_score(y_test, y_scores)
    baseline = np.bincount(y_test)[1] / sum(np.bincount(y_test))
    plt.figure()
    plt.step(recall, precision, color='b', alpha=0.2,
             where='post')
    plt.fill_between(recall, precision, step='post', alpha=0.2,
                     color='b')
    plt.axhline(y=baseline, linewidth=2, color='navy', linestyle='--')
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.ylim([0.0, 1.05])
    plt.xlim([0.0, 1.0])
    plt.title('2-class Precision-Recall curve: AP={0:0.2f}'.format(
              average_precision))
    plt.show()


def reportclassfication(clf, X_test, y_test):
    '''Reports classification results with the given model and testdata'''

    print("Detailed classification report:")
    print()
    y_true, y_pred = y_test, clf.predict(X_test)
    print(classification_report(y_true, y_pred))
    print()


def timeit(method):
    """Decorator for time measurement"""
    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()
        print ('%r  %2.2f ms' % \
            (method.__name__, (te - ts) * 1000))
        return result
    return timed


@timeit
def apply_motif_kernel(motifs, sequences,return_kernel_matrix = False):
    return motifKernel(motifs, sequences, return_kernel_matrix)


def main():
    # load the data
    # stemcells
    pos_data = [seq.seq for seq in sio.parse('data/stemcells.fa', 'fasta')]
    # fibroblasts
    neg_data = [seq.seq for seq in sio.parse('data/fibroblast.fa', 'fasta')]

    pos_labels = np.ones(len(pos_data), dtype=int)
    neg_labels = np.zeros(len(neg_data), dtype=int)

    y = np.concatenate((pos_labels, neg_labels), axis=0)

    with open("sequences.txt","w") as testseqs:
        for seq in pos_data + neg_data:
            testseqs.write(str(seq) + "\n")

    # motifs
    motifs = ["TCAGCA","TGCTGA","ATGCA.A","T.TGCAT"]
    print("Number of Motifs: {}".format(len(motifs)))
    # motif kernel and print execution time
    X = apply_motif_kernel(motifs, pos_data + neg_data,False)

    # split data
    X_train, X_test, y_train, y_test = train_test_split(X,
                                                        y,
                                                        test_size=0.1,
                                                        random_state=42,
                                                        stratify=y)
    #train classifier
    clf = SVC()
    clf.fit(X_train, y_train)

    # get classfication results with the tuned parameters
    reportclassfication(clf, X_test, y_test)

    y_scores = clf.decision_function(X_test)
    # plot ROC
    plot_roc_curve(y_test, y_scores)
    # plot PRC
    plot_prec_recall_curve(y_test, y_scores)


if __name__ == "__main__":
    main()
