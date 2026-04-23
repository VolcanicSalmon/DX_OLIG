for i in *hopag1_*/*cif  
do              
python3 ../multidom_contacts.py -c $i -d 'VIETDKWSQLSGAKGSNPGGLFQAPDGVKWYVKTNPSTNRLRNEVLASKLYRAAGIDVPDIKMASRQGKPALISKLIGGNHKDIKTIEGSGQLRSGFAVDAWLANWDVIGQTGDNIIFNDRNKPVRIDLGGSLLFRAQGGHKGNQFGNTPMELVTMLSRKENTSSHAFRKIERNDIRMGIAAIERIPDDRIKALCAEHGPGNYSERIELGKRLISRKHWLVNMKQTLPHIHRQK' -o $(basename $i)_Kinase
done
for i in ar6dock2/at2g26400.2_model ar6dock2/at2g26400.3_model ar6dock2/at5g48410.1_model ar6dock2/at5g48410.2_model ar6dock/at2g26400.4_model ar6dock/at4g22880.3_model ar6dock/at5g44990.2_model; do     ^Cthon3 hadconf_prot*2.py -r $i -p prots/wt -t /hpc-home/vef25hok/DANXIA/ligands4dock/AR6_ideal.acpype/AR6_ideal_CNS.top \
-P /hpc-home/vef25hok/DANXIA/ligands4dock/AR6_ideal.acpype/AR6_ideal_CNS.par -o wt_ar6; done
