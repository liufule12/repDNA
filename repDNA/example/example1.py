__author__ = 'liufule12'

import time

from repDNA.nac import RevcKmer
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

    # Generate the feature vectors based on reverse compliment kmer.
    rev_kmer = RevcKmer(k=6, normalize=True, upto=True)
    pos_vec = rev_kmer.make_revckmer_vec(open('hs.fasta'))
    neg_vec = rev_kmer.make_revckmer_vec(open('non-hs.fasta'))

    print(len(pos_vec))
    print(len(neg_vec))

    # Merge positive and negative feature vectors and generate their corresponding labels.
    vec = np.array(pos_vec + neg_vec)
    vec_label = np.array([0] * len(pos_vec) + [1] * len(neg_vec))

    # ##############################################################################
    # Classification and accurate analysis.

    # Using 10-fold cross-validation to evaluate the performance of the predictor.
    clf = svm.LinearSVC()
    scores = cross_validation.cross_val_score(clf, vec, y=vec_label, cv=10)
    print('Per accuracy in 10-fold CV:')
    print(scores)
    print("Accuracy: %0.2f (+/- %0.2f)" % (scores.mean(), scores.std() * 2))

    # ###############################################################################
    # Classification and ROC analysis.

    # Run classifier with cross-validation and plot ROC curves
    cv = StratifiedKFold(vec_label, n_folds=10)
    classifier = svm.SVC(kernel='linear', probability=True)

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

    print('Example1 End.')

    total_time = time.time() - begin_time
    print('Total running time of the example: %.2f seconds ( %i minutes %.2f seconds )' % (
        total_time, int(total_time / 60), total_time % 60))