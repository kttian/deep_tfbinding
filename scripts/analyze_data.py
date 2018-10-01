from __future__ import print_function

DataDir = "./"
header_file = DataDir + "headers.txt"
train_file  = DataDir + "finetune/splits/train.tsv.gz"
valid_file  = DataDir + "splits/valid.tsv.gz"
test_file   = DataDir + "splits/test.tsv.gz"

import gzip
import pandas as pd
import numpy as np
from tabulate import tabulate

def analyze_file(in_name, out_name, head):

    df0 = pd.read_csv(in_name, sep='\t', index_col=None, header=None, compression='gzip')
    df = df0.iloc[:,3:]
    df.columns = head

    counts_list = []
    for c in range(df.shape[1]):
        counts = df.iloc[:, c].value_counts().sort_index(ascending=False)
        counts_list.append(counts)

    counts_df = pd.concat(counts_list, axis=1).sort_index(ascending=False)
    #print(counts_df)
    sums = counts_df.sum()
    #print (sums)

    percent_df = counts_df / sums

    new_index = [ int(1), int(0), int(-1)]
    counts_df  = counts_df.reindex(new_index)
    percent_df = percent_df.reindex(new_index)

    all=pd.concat([counts_df, percent_df])

    all_t = all.transpose()
    all_t.columns=['pos', 'neg', 'amb', 'pos%', 'neg%', 'amb%']
    all_t.to_csv(out_name, sep='\t')

    print(tabulate(all_t, headers='keys', tablefmt='psql'))

def draw_roc_curves(prefix, all_y, all_predictions, task_id=-1, tf_name=""):
    import matplotlib as mpl
    mpl.use("Agg")
    import matplotlib.pyplot as plt
    from sklearn.metrics import roc_curve, auc

    n_classes = all_y.shape[1]

    # Compute ROC curve and ROC area for each class
    fpr = dict()
    tpr = dict()
    roc_auc = dict()
    for i in range(n_classes):
        if task_id != -1 and i != task_id:
            continue
        true_y_for_task = all_y[:,i]
        predictions_for_task = all_predictions[:,i]
        predictions_for_task_filtered, true_y_for_task_filtered = \
         remove_ambiguous_peaks(predictions_for_task, true_y_for_task)
        fpr[i], tpr[i], _ = roc_curve(true_y_for_task_filtered,
                                      predictions_for_task_filtered)
        roc_auc[i] = auc(fpr[i], tpr[i])
    plt.figure()
    lw=1
    color_list = np.linspace(0, 1, n_classes)

    for i in range(n_classes):
        color = plt.cm.plasma(color_list[i])
        if task_id != -1 and i != task_id:
            continue
        plt.plot(fpr[i], tpr[i], color=color, lw=lw,
                 label='ROC curve of class {0} (area = {1:0.2f})'
                 ''.format(i, roc_auc[i]))
    plt.plot([0, 1], [0, 1], 'k--', lw=lw)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(tf_name + ' Receiver Operating Characteristics')
    plt.legend(loc="lower right")
    #plt.show()
    plt.savefig(prefix + "ROC.png")

def draw_prc_curves(prefix, all_y, all_predictions, task_id=-1, tf_name=""):
    import matplotlib as mpl
    mpl.use("Agg")
    import matplotlib.pyplot as plt
    from sklearn.metrics import precision_recall_curve, auc, average_precision_score

    n_classes = all_y.shape[1]

    # Compute ROC curve and ROC area for each class
    precision = dict()
    recall = dict()
    prc_auc = dict()
    for i in range(n_classes):
        if task_id != -1 and i != task_id:
            continue
        true_y_for_task = all_y[:,i]
        predictions_for_task = all_predictions[:,i]
        predictions_for_task_filtered, true_y_for_task_filtered = \
         remove_ambiguous_peaks(predictions_for_task, true_y_for_task)
        precision[i], recall[i], _ = precision_recall_curve(true_y_for_task_filtered,
                                                            predictions_for_task_filtered)
        prc_auc[i] = auc(recall[i], precision[i])
        prc_auc[i] = average_precision_score(true_y_for_task_filtered, predictions_for_task_filtered)

    plt.figure()
    lw=1
    color_list = np.linspace(0, 1, n_classes)

    for i in range(n_classes):
        color = plt.cm.plasma(color_list[i])
        if task_id != -1 and i != task_id:
            continue
        plt.plot(recall[i], precision[i], color=color, lw=lw,
                 label='PRC curve of class {0} (area = {1:0.2f})'
                 ''.format(i, prc_auc[i]))
    plt.plot([0, 1], [1, 0], 'k--', lw=lw)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title(tf_name + ' Precision-recall Curve')
    plt.legend(loc="lower left")
    #plt.show()
    plt.savefig(prefix + "PRC.png")


def load_predictions(predict_file):
    pred = np.load(predict_file)
    assert(pred.shape[1] % 2 ==0)
    n_cols = int(pred.shape[1] / 2)
    all_y = pred[:, :n_cols]
    all_predictions = pred[:, n_cols:]
    return all_y, all_predictions

# Python script for confusion matrix creation.
from pandas_ml import ConfusionMatrix

from sklearn.metrics import confusion_matrix
from sklearn.metrics import accuracy_score
from sklearn.metrics import classification_report

def remove_ambiguous_peaks(true_y, predictions): 
    indices_to_remove = np.where(true_y < 0)
    true_y_filtered = np.delete(true_y, indices_to_remove)
    predictions_filtered = np.delete(predictions, indices_to_remove)
    return true_y_filtered, predictions_filtered

def analyze_result(in_y, in_predictions):
    
    y, predict_prob = remove_ambiguous_peaks(in_y, in_predictions)

    predictions = np.where(predict_prob > 0.5, 1, 0)

    #cm = confusion_matrix(y, predictions)
    #print('Confusion Matrix :')
    #print(cm)
    cm = ConfusionMatrix(y, predictions)
    print("Confusion matrix:\n%s\n" % cm)
    print('Accuracy Score :',accuracy_score(y, predictions))
    print('Report : ')
    print(classification_report(y, predictions))

with open(header_file, 'r') as f:
    head = f.readline()
    head = head.split()[1:]

print("valid set")
analyze_file(valid_file, "valid_counts.tsv", head)
print("test set")
analyze_file(test_file,  "test_counts.tsv", head)
print("train set")
analyze_file(train_file, "train_counts.tsv", head)

all_y, all_predictions = load_predictions("model_files/record_1_predictions.npy")

num_tasks = all_y.shape[1]

for t in range(num_tasks):
    print("")
    print("task %d %s -------" % (t, head[t]))
    analyze_result(all_y[:, t], all_predictions[:, t])
