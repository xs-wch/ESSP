function auc = myauc(fpr,tpr)
auc = sum(diff([0,fpr]).*tpr);
end