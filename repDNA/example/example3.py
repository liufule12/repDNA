__author__ = 'Fule Liu'

import time

from repDNA.psenac import PseDNC
import numpy as np
from sklearn import svm
from sklearn import cross_validation
from sklearn.cross_validation import StratifiedKFold
from sklearn.metrics import roc_curve, auc
from scipy import interp
import matplotlib.pyplot as plt


if __name__ == '__main__':
    begin_time = time.time()
    print('Example1 Start.(This process may use several minutes, please do not close the program.)')

    # ##############################################################################
    # Data IO and generation.

    # Generate the PseDNC feature vector.
    psednc = PseDNC(lamada=3, w=0.05)
    pos_vec = psednc.make_psednc_vec(open('hotspots.fasta'))
    neg_vec = psednc.make_psednc_vec(open('coldspots.fasta'))

    print(len(pos_vec))
    print(len(neg_vec))

    # Merge positive and negative feature vectors and generate their corresponding labels.
    vec = np.array(pos_vec + neg_vec)
    vec_label = np.array([1] * len(pos_vec) + [0] * len(neg_vec))

    # ##############################################################################
    # Classification and accurate analysis.

    # evaluate performance of the predictor by 5-fold cross-validation and plot the mean ROC curve.
    clf = svm.SVC(C=32, gamma=0.5)
    scores = cross_validation.cross_val_score(clf, vec, y=vec_label, cv=5)
    print('Per accuracy in 5-fold CV:')
    print(scores)
    print("Accuracy: %0.2f (+/- %0.2f)" % (scores.mean(), scores.std() * 2))

    # ###############################################################################
    # Classification and ROC analysis.

    # evaluate performance of the predictor by 5-fold cross-validation and plot the mean ROC curve.
    cv = StratifiedKFold(vec_label, n_folds=5)
    classifier = svm.SVC(C=32, kernel='rbf', gamma=0.5, probability=True)

    mean_tpr = 0.0
    mean_fpr = np.linspace(0, 1, 100)
    all_tpr = []

    for i, (train, test) in enumerate(cv):
        probas_ = classifier.fit(vec[train], vec_label[train]).predict_proba(vec[test])
        # Compute ROC curve and area the curve
        fpr, tpr, thresholds = roc_curve(vec_label[test], probas_[:, 1])
        mean_tpr += interp(mean_fpr, fpr, tpr)
        mean_tpr[0] = 0.0

    # Plot ROC curve.
    mean_tpr /= len(cv)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    plt.plot(mean_fpr, mean_tpr, '-', label='Mean ROC (area = %0.2f)' % mean_auc, lw=2)

    plt.xlim([0, 1.0])
    plt.ylim([0, 1.0])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.legend(loc="lower right")
    plt.show()

    print('Example3 End.')

    total_time = time.time() - begin_time
    print('Total running time of the example: %.2f seconds ( %i minutes %.2f seconds )' % (
        total_time, int(total_time / 60), total_time % 60))