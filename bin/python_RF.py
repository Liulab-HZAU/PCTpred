#!/usr/bin/env python
from sklearn.ensemble import RandomForestClassifier
import numpy, argparse, sys

#==============================Classification trees and forests======================

def read_samples(file):
    feature_data = []
    labels = []
    for line in open(file):
        line = line.rstrip()
        records = line.split('\t')
        labels.append(records[0])
        values = [float(x) for x in records[1:]]
        feature_data.append(values)
    return feature_data,labels

def main(args = sys.argv):
    parser = argparse.ArgumentParser('Random Forest Classifier and tree interpreter')
    parser.add_argument('-train', '--train', default='train.txt',
                        help='train file')
    parser.add_argument('-test', '--test', default='test.txt',
                        help='test file')
    parser.add_argument('-o', '--out', default='out.txt',
                        help='out file')
    args = parser.parse_args()
    
    train_features,train_labels = read_samples(args.train)
    test_features,test_labels = read_samples(args.test)
    
    rf = RandomForestClassifier(n_estimators= 500)
    rf.fit(train_features,train_labels)
    prediction = rf.predict_proba(test_features)

    f = open(args.out, 'w')
    f.write('0\t1\n')
    n = 0 
    for i in prediction:
        n = n + 1
        f.write('{0}\t{1}\t{2}\n'.format(n,i[0],i[1]))
    f.close()
if __name__ == '__main__': 
    main()