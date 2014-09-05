__author__ = 'liufule12'

import time

from repDNA.ac import DAC
import numpy as np
from sklearn import svm
from sklearn import cross_validation
from sklearn.cross_validation import StratifiedKFold
from sklearn.metrics import roc_curve, auc
from scipy import interp
import matplotlib.pyplot as plt


if __name__ == '__main__':
    begin_time = time.time()
    print 'Example1 Start.(This process may use several minutes, please do not close the program.)'

    # ##############################################################################
    # Data IO and generation.

    # Generate the corresponding feature vector.
    ac = DAC(lag=6)
    pos_vec = ac.make_dac_vec(open('H_sapiens_pos.fasta'), all_property=True)
    neg_vec = ac.make_dac_vec(open('H_sapiens_neg.fasta'), all_property=True)

    print 'The number of positive and negative vectors.'
    print len(pos_vec)
    print len(neg_vec)

    print 'The dimension of per positive and negative vector'
    print len(pos_vec[0])
    print len(neg_vec[0])

    # write_libsvm(pos_vec, neg_vec, write_filename)

    # Merge feature vector and generate corresponding vector label.
    vec = np.array(pos_vec + neg_vec)
    vec_label = np.array([0] * len(pos_vec) + [1] * len(neg_vec))

    # ##############################################################################
    # Classification and accurate analysis.

    # Run classifier with 10 folds cross-validation and generate accuracy.
    clf = svm.SVC(C=32768.0)
    scores = cross_validation.cross_val_score(clf, vec, y=vec_label, cv=5)
    print 'Per accuracy in 5-fold CV:'
    print scores
    print("Accuracy: %0.2f (+/- %0.2f)" % (scores.mean(), scores.std() * 2))

    # ###############################################################################
    # Classification and ROC analysis.

    # Run classifier with cross-validation and plot ROC curves
    cv = StratifiedKFold(vec_label, n_folds=5)
    classifier = svm.SVC(C=32768.0, kernel='rbf', gamma=0.001953125,
                         probability=True)

    mean_tpr = 0.0
    mean_fpr = np.linspace(0, 1, 100)
    all_tpr = []

    for i, (train, test) in enumerate(cv):
        probas_ = classifier.fit(vec[train],
                                 vec_label[train]).predict_proba(vec[test])
        # Compute ROC curve and area the curve
        fpr, tpr, thresholds = roc_curve(vec_label[test], probas_[:, 1])
        mean_tpr += interp(mean_fpr, fpr, tpr)
        mean_tpr[0] = 0.0

    # Plot ROC curve.
    plt.plot([0, 1], [0, 1], '--', color=(0.6, 0.6, 0.6), label='Luck')

    mean_tpr /= len(cv)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    plt.plot(mean_fpr, mean_tpr, '-',
             label='Mean ROC (area = %0.2f)' % mean_auc, lw=2)

    plt.xlim([0, 1.0])
    plt.ylim([0, 1.0])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    # plt.title('ROC curve for the 5-fold cross-validation of DAC')
    plt.legend(loc="lower right")
    plt.show()

    print 'Example2 End.'

    total_time = time.time() - begin_time
    print 'Total running time of the example: %.2f seconds ( %i minutes %.2f seconds )' % (
        total_time, int(total_time / 60), total_time % 60)