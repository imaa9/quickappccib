#!/bin/bash

#SBATCH -p shared
#SBATCH --cpus-per-task=4
#SBATCH -N 1
#SBATCH --mem 16G
#SBATCH -t 3-00:00:00
#SBATCH -J PCAs
#SBATCH -o PCAs_%j.out
#SBATCH -e PCAs_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=imaayan@g.harvard.edu
#SBATCH --export=NONE

module load Mambaforge/22.11.1-fasrc01
mamba activate ipyrad

cd /n/holyscratch01/losos_lab/Users/imaayan/all_ddrad/Analysis/PCA

# grahami 
python -c'

import ipyrad.analysis as ipa
import pandas as pd
import toyplot

data = "/n/holyscratch01/losos_lab/Users/imaayan/all_ddrad/ipyrad/grahami/grh_G95_outfiles/grh_G95.snps.hdf5"

imap={
"grahami":["grh_14581","grh_14591","grh_14578","grh_149864","grh_149860","grh_149958","grh_14634","grh_14655","grh_14656","grh_150384","grh_150390","grh_150714","grh_150715","grh_150613","grh_150608","grh_150615","grh_198955","grh_150057","grh_150069","grh_150060","grh_150059","grh_150077","grh_150079","grh_150086","grh_150083","grh_150799","grh_150798","grh_149894","grh_149892","grh_14614","grh_14605","grh_14602","grh_150458","grh_2674","grh_2641","grh_14784","grh_14788","grh_198958","grh_150524","grh_150521","grh_149849","grh_2581","grh_2659","grh_2658","grh_2637","grh_2663","grh_2639","grh_14710","grh_14723","grh_14722","grh_14706","grh_14715","grh_149950","grh_149943","grh_150625","grh_150650","grh_3707","grh_3683","grh_3882ra","grh_150765","grh_150188","grh_150192","grh_150416","grh_149909","grh_150170","grh_150169","grh_150282","grh_150290","grh_150814","grh_150482","grh_150174","grh_150179","grh_150606","grh_150037","grh_2613","grh_150260","grh_150245","grh_150154","grh_150141","grh_150783","grh_14552","grh_14727","grh_150264","grh_150261","grh_14740","grh_14741","grh_14742","grh_14744","grh_2606","grh_150441","grh_150534","grh_150537","grh_14560","grh_14561","grh_14570","grh_14559","grh_14557","grh_149818","grh_149825","grh_150496","grh_150497","grh_150377","grh_14677","grh_14658","grh_150322","grh_150301","grh_14680","grh_14679","grh_150669","grh_150668","grh_2760","grh_2747","grh_14800","grh_14794","grh_150008","grh_14808","grh_14797","grh_14796","grh_150116","grh_150109","grh_150104","grh_150101","grh_150205","grh_150203","grh_149967","grh_150098","grh_150587","grh_150583","grh_150581","grh_150586","grh_198962","grh_149868","grh_149866","grh_149869","grh_150200","grh_150193","grh_150223","grh_150229","grh_150234","grh_2622","grh_2577","grh_2624","grh_2623","grh_150746","grh_150758","grh_150753","grh_150699","grh_150694","grh_150343","grh_150340","grh_150810","grh_150808","grh_14826","grh_14824","grh_14825","grh_150427","grh_149899","grh_149896","grh_149898","grh_150687","grh_150596","grh_150591","grh_150604","grh_150168","grh_150162","grh_150163","grh_14537","grh_14536","grh_14539","grh_150575","grh_150564"],
}

pca = ipa.pca(
data=data,
imap=imap,
minmap=0,
mincov=0.95,
impute_method="sample",
)

pca.run(nreplicates=100, seed=2140)
df = pd.DataFrame(pca.pcaxes[0], index=pca.names)
df.to_csv("pca_grahami95.csv")
dfv = pd.DataFrame(pca.variances)
dfv.to_csv("var_grahami95.csv")
dfm = pd.DataFrame(pca.missing)
dfm.to_csv("miss_grahami95.csv")
'

sed -i "1s/,/,PC/g" pca_grahami95.csv
sed -i "1s/^/ID/g" pca_grahami95.csv
sed -i "1s/,/,rep/g" var_grahami95.csv
sed -i "1s/^/PC/g" var_grahami95.csv
sed -i "1s/^/ID/g" miss_grahami95.csv


# opalinus 
python -c'

import ipyrad.analysis as ipa
import pandas as pd
import toyplot

data = "/n/holyscratch01/losos_lab/Users/imaayan/all_ddrad/ipyrad/opalinus/opl_G95_outfiles/opl_G95.snps.hdf5"

imap={
"opalinus":["opl_14587", "opl_14584", "opl_14583", "opl_14580ra", "opl_14577", "opl_149861", "opl_149960", "opl_149957", "opl_14642", "opl_14633", "opl_14650", "opl_14652", "opl_150055", "opl_150056", "opl_150397", "opl_150398", "opl_14613", "opl_150454", "opl_150456", "opl_14781", "opl_14780", "opl_14785", "opl_149934", "opl_149938", "opl_150520", "opl_150523", "opl_150525", "opl_198965", "opl_150532", "opl_150513", "opl_150511", "opl_14708", "opl_14714", "opl_14704", "opl_14720", "opl_149949", "opl_149944", "opl_150766", "opl_150763", "opl_150403", "opl_150408", "opl_149978", "opl_149979", "opl_150481", "opl_150484", "opl_150607", "opl_150034", "opl_150371", "opl_150367", "opl_150242", "opl_150243", "opl_150774", "opl_150788", "opl_14554", "opl_150218", "opl_150220", "opl_150217", "opl_198963", "opl_14752", "opl_14764", "opl_14756", "opl_14758", "opl_150215", "opl_2676", "opl_150538", "opl_150553", "opl_150378", "opl_14617", "opl_14621", "opl_14624", "opl_14626", "opl_150312", "opl_150302", "opl_14798", "opl_150012", "opl_150009", "opl_14789", "opl_14795", "opl_14791", "opl_14792", "opl_150113", "opl_150103", "opl_150202", "opl_149966", "opl_149965", "opl_149996", "opl_150007", "opl_150093", "opl_150092", "opl_149867", "opl_198960", "opl_150041", "opl_150045", "opl_150044", "opl_150046", "opl_14664", "opl_14667", "opl_14666", "opl_150421", "opl_14775", "opl_14776", "opl_14778", "opl_14777", "opl_149983", "opl_149981", "opl_150738", "opl_150748", "opl_150809", "opl_14827", "opl_14830", "opl_14828", "opl_150470", "opl_150467", "opl_150131", "opl_150130", "opl_150128", "opl_14837", "opl_14838", "opl_14814", "opl_150125", "opl_150123", "opl_150121", "opl_14533", "opl_14535", "opl_14529", "opl_14528", "opl_14526", "opl_150558", "opl_150573","opl_149908","opl_149916","opl_149927","opl_149928",], 
}

pca = ipa.pca(
data=data,
imap=imap,
minmap=0,
mincov=0.95,
impute_method="sample",
)

pca.run(nreplicates=100, seed=2141)
df = pd.DataFrame(pca.pcaxes[0], index=pca.names)
df.to_csv("pca_opalinus95.csv")
dfv = pd.DataFrame(pca.variances)
dfv.to_csv("var_opalinus95.csv")
dfm = pd.DataFrame(pca.missing)
dfm.to_csv("miss_opalinus95.csv")
'

sed -i "1s/,/,PC/g" pca_opalinus95.csv
sed -i "1s/^/ID/g" pca_opalinus95.csv
sed -i "1s/,/,rep/g" var_opalinus95.csv
sed -i "1s/^/PC/g" var_opalinus95.csv
sed -i "1s/^/ID/g" miss_opalinus95.csv



# lineatopus 
python -c'

import ipyrad.analysis as ipa
import pandas as pd
import toyplot

data = "/n/holyscratch01/losos_lab/Users/imaayan/all_ddrad/ipyrad/lineatopus/ltp_G95_outfiles/ltp_G95.snps.hdf5"

imap={
"lineatopus":["ltp_14589", "ltp_14592", "ltp_14576", "ltp_149863", "ltp_149859", "ltp_150331", "ltp_150330", "ltp_150329", "ltp_150353", "ltp_150361", "ltp_14657", "ltp_14654", "ltp_150382", "ltp_150387", "ltp_150381", "ltp_150716", "ltp_150712", "ltp_150619", "ltp_150617", "ltp_150610", "ltp_150064", "ltp_150068", "ltp_150070", "ltp_150072", "ltp_150393", "ltp_149888", "ltp_149893", "ltp_14603", "ltp_14600", "ltp_14608", "ltp_14607", "ltp_14782", "ltp_14786", "ltp_14787", "ltp_150527", "ltp_150529", "ltp_150530", "ltp_150514", "ltp_149850", "ltp_149851", "ltp_149856", "ltp_149857", "ltp_149874", "ltp_149878", "ltp_149872", "ltp_14709", "ltp_14705", "ltp_149947", "ltp_149946", "ltp_150682", "ltp_150628", "ltp_150643", "ltp_150640", "ltp_3829", "ltp_3469", "ltp_3494", "ltp_150767", "ltp_198964", "ltp_150187", "ltp_150189", "ltp_150407", "ltp_150413", "ltp_149920", "ltp_149914", "ltp_149907", "ltp_149977", "ltp_150284", "ltp_150289", "ltp_150479", "ltp_150487", "ltp_149929", "ltp_149930", "ltp_149925", "ltp_149931", "ltp_14692", "ltp_14688", "ltp_14684", "ltp_14683", "ltp_14682rb", "ltp_150175", "ltp_150172", "ltp_150032", "ltp_150030", "ltp_150374", "ltp_150248", "ltp_150240", "ltp_150258", "ltp_150152", "ltp_150150", "ltp_150780", "ltp_150797", "ltp_14548", "ltp_14543", "ltp_14544", "ltp_14550", "ltp_14540", "ltp_14731", "ltp_14733", "ltp_150263", "ltp_150265", "ltp_14736", "ltp_14739", "ltp_150436", "ltp_150431", "ltp_150539", "ltp_14572", "ltp_14558", "ltp_14566", "ltp_14562", "ltp_198966", "ltp_149822", "ltp_149820", "ltp_149831", "ltp_150501", "ltp_150499", "ltp_149989", "ltp_150375", "ltp_14629", "ltp_14628", "ltp_150324", "ltp_150327", "ltp_150305", "ltp_150325", "ltp_150303", "ltp_150304", "ltp_14676", "ltp_150666", "ltp_150670", "ltp_14804", "ltp_14801", "ltp_14805", "ltp_150018", "ltp_150019", "ltp_14810", "ltp_14809", "ltp_150110", "ltp_150115", "ltp_150111", "ltp_150100", "ltp_150105", "ltp_150106", "ltp_150212", "ltp_150206", "ltp_149961", "ltp_149993", "ltp_150006", "ltp_150094", "ltp_150090", "ltp_149971", "ltp_149969", "ltp_149972", "ltp_149975", "ltp_150049", "ltp_150052", "ltp_150042", "ltp_150198", "ltp_150196", "ltp_150239", "ltp_150221", "ltp_150222", "ltp_150238", "ltp_14670", "ltp_14671", "ltp_14663", "ltp_14668", "ltp_14665", "ltp_14669", "ltp_150425", "ltp_150426", "ltp_150417", "ltp_150428", "ltp_14773", "ltp_198957", "ltp_149982", "ltp_149985", "ltp_149988", "ltp_14765", "ltp_14772", "ltp_14768", "ltp_14769", "ltp_150736", "ltp_150739", "ltp_150751", "ltp_150754", "ltp_150728", "ltp_150724", "ltp_150704", "ltp_150696", "ltp_149842", "ltp_149840", "ltp_150339", "ltp_150342", "ltp_150806", "ltp_150812", "ltp_14822", "ltp_14823", "ltp_14821", "ltp_14820", "ltp_150707", "ltp_150709", "ltp_149882", "ltp_149880", "ltp_149883", "ltp_149901", "ltp_149904", "ltp_149902", "ltp_149903", "ltp_150469", "ltp_150465", "ltp_150462", "ltp_150685", "ltp_149953", "ltp_149955", "ltp_149954", "ltp_149956", "ltp_150602", "ltp_150600", "ltp_14836", "ltp_14834", "ltp_14833", "ltp_150117", "ltp_14815", "ltp_14817", "ltp_150122", "ltp_150118", "ltp_150160", "ltp_14530", "ltp_14534", "ltp_14522", "ltp_14523", "ltp_198961", "ltp_150570", "ltp_150566", "ltp_150560"],
}

pca = ipa.pca(
data=data,
imap=imap,
minmap=0,
mincov=0.95,
impute_method="sample",
)

pca.run(nreplicates=100, seed=2142)
df = pd.DataFrame(pca.pcaxes[0], index=pca.names)
df.to_csv("pca_lineatopus95.csv")
dfv = pd.DataFrame(pca.variances)
dfv.to_csv("var_lineatopus95.csv")
dfm = pd.DataFrame(pca.missing)
dfm.to_csv("miss_lineatopus95.csv")
'

sed -i "1s/,/,PC/g" pca_lineatopus95.csv
sed -i "1s/^/ID/g" pca_lineatopus95.csv
sed -i "1s/,/,rep/g" var_lineatopus95.csv
sed -i "1s/^/PC/g" var_lineatopus95.csv
sed -i "1s/^/ID/g" miss_lineatopus95.csv


# valencienni 
python -c'

import ipyrad.analysis as ipa
import pandas as pd
import toyplot

data = "/n/holyscratch01/losos_lab/Users/imaayan/all_ddrad/ipyrad/valencienni/vnc_G95_outfiles/vnc_G95.snps.hdf5"

imap={
"valencienni":["vnc_150076", "vnc_150078", "vnc_150081", "vnc_150653", "vnc_150629", "vnc_150642", "vnc_3546", "vnc_3509", "vnc_3533", "vnc_3559", "vnc_3573", "vnc_150191", "vnc_150412", "vnc_14698", "vnc_14689", "vnc_14685", "vnc_14699", "vnc_150368", "vnc_150369", "vnc_150241", "vnc_150793", "vnc_14738", "vnc_14737", "vnc_150271", "vnc_150273", "vnc_14743", "vnc_150554", "vnc_198953", "vnc_150556", "vnc_14573", "vnc_14632", "vnc_14674", "vnc_14678", "vnc_14681", "vnc_14660", "vnc_14661", "vnc_14675", "vnc_150011", "vnc_14812", "vnc_14802", "vnc_150102", "vnc_150000", "vnc_150227", "vnc_150232", "vnc_150228", "vnc_14662", "vnc_149834", "vnc_150695", "vnc_14829", "vnc_14831"],
}

pca = ipa.pca(
data=data,
imap=imap,
minmap=0,
mincov=0.95,
impute_method="sample",
)

pca.run(nreplicates=100, seed=2143)
df = pd.DataFrame(pca.pcaxes[0], index=pca.names)
df.to_csv("pca_valencienni95.csv")
dfv = pd.DataFrame(pca.variances)
dfv.to_csv("var_valencienni95.csv")
dfm = pd.DataFrame(pca.missing)
dfm.to_csv("miss_valencienni95.csv")
'

sed -i "1s/,/,PC/g" pca_valencienni95.csv
sed -i "1s/^/ID/g" pca_valencienni95.csv
sed -i "1s/,/,rep/g" var_valencienni95.csv
sed -i "1s/^/PC/g" var_valencienni95.csv
sed -i "1s/^/ID/g" miss_valencienni95.csv


# conspersus 
python -c'

import ipyrad.analysis as ipa
import pandas as pd
import toyplot

data = "/n/holyscratch01/losos_lab/Users/imaayan/all_ddrad/ipyrad/conspersus/cnp_G95_outfiles/cnp_G95.snps.hdf5"

imap={
"conspersus":["cnp_8253", "cnp_8249", "cnp_8254", "cnp_8279", "cnp_8280", "cnp_8283", "cnp_8218ra", "cnp_8031", "cnp_8270", "cnp_8227", "cnp_8139", "cnp_8145", "cnp_8159", "cnp_8208", "cnp_8211", "cnp_8213", "cnp_8240", "cnp_8239", "cnp_8242", "cnp_8261", "cnp_8306", "cnp_8292"],
}

pca = ipa.pca(
data=data,
imap=imap,
minmap=0,
mincov=0.95,
impute_method="sample",
)

pca.run(nreplicates=100, seed=2144)
df = pd.DataFrame(pca.pcaxes[0], index=pca.names)
df.to_csv("pca_conspersus95.csv")
dfv = pd.DataFrame(pca.variances)
dfv.to_csv("var_conspersus95.csv")
dfm = pd.DataFrame(pca.missing)
dfm.to_csv("miss_conspersus95.csv")
'

sed -i "1s/,/,PC/g" pca_conspersus95.csv
sed -i "1s/^/ID/g" pca_conspersus95.csv
sed -i "1s/,/,rep/g" var_conspersus95.csv
sed -i "1s/^/PC/g" var_conspersus95.csv
sed -i "1s/^/ID/g" miss_conspersus95.csv


# garmani 
python -c'

import ipyrad.analysis as ipa
import pandas as pd
import toyplot

data = "/n/holyscratch01/losos_lab/Users/imaayan/all_ddrad/ipyrad/garmani/gmn_G95_outfiles/gmn_G95.snps.hdf5"

imap={
"garmani":["gmn_14593", "gmn_14588", "gmn_14594", "gmn_14582", "gmn_150352", "gmn_150354", "gmn_14638", "gmn_14643", "gmn_14641", "gmn_14646", "gmn_150392", "gmn_150380", "gmn_150711", "gmn_150720", "gmn_150623", "gmn_150054", "gmn_150089", "gmn_150087", "gmn_150400", "gmn_149886", "gmn_14606", "gmn_14597", "gmn_14599", "gmn_14609", "gmn_14604", "gmn_150451", "gmn_150449", "gmn_149937",  "gmn_149843","gmn_149858", "gmn_149875", "gmn_149876", "gmn_14711", "gmn_14718ra", "gmn_150678", "gmn_150684", "gmn_149952", "gmn_150655", "gmn_150626", "gmn_150761", "gmn_150759", "gmn_149905", "gmn_198954", "gmn_150294", "gmn_149923", "gmn_14686", "gmn_14696", "gmn_14694", "gmn_14693", "gmn_150605", "gmn_150036", "gmn_150035", "gmn_150365", "gmn_150795", "gmn_150794", "gmn_14547", "gmn_14551", "gmn_14763", "gmn_14757", "gmn_14755", "gmn_150536", "gmn_150555", "gmn_14568", "gmn_14563", "gmn_14569", "gmn_149990", "gmn_14618", "gmn_14623", "gmn_150315", "gmn_150311", "gmn_150673", "gmn_149951", "gmn_150112", "gmn_150114", "gmn_149962", "gmn_149998", "gmn_149991", "gmn_150585", "gmn_150578", "gmn_14673", "gmn_14672", "gmn_150735", "gmn_150752", "gmn_150725", "gmn_150697", "gmn_149846", "gmn_150811", "gmn_150807", "gmn_150708", "gmn_149881", "gmn_150472", "gmn_150468", "gmn_150595", "gmn_150599", "gmn_14832", "gmn_150129", "gmn_150126", "gmn_14819", "gmn_14816rb", "gmn_14818", "gmn_150120", "gmn_150166"], 
}

pca = ipa.pca(
data=data,
imap=imap,
minmap=0,
mincov=0.95,
impute_method="sample",
)

pca.run(nreplicates=100, seed=2145)
df = pd.DataFrame(pca.pcaxes[0], index=pca.names)
df.to_csv("pca_garmani95.csv")
dfv = pd.DataFrame(pca.variances)
dfv.to_csv("var_garmani95.csv")
dfm = pd.DataFrame(pca.missing)
dfm.to_csv("miss_garmani95.csv")
'

sed -i "1s/,/,PC/g" pca_garmani95.csv
sed -i "1s/^/ID/g" pca_garmani95.csv
sed -i "1s/,/,rep/g" var_garmani95.csv
sed -i "1s/^/PC/g" var_garmani95.csv
sed -i "1s/^/ID/g" miss_garmani95.csv

# reconditus 
python -c'

import ipyrad.analysis as ipa
import pandas as pd
import toyplot

data = "/n/holyscratch01/losos_lab/Users/imaayan/all_ddrad/ipyrad/reconditus/rcn_G95_outfiles/rcn_G95.snps.hdf5"

imap={"reconditus":["rcn_150133", "rcn_150136", "rcn_198956", "rcn_150139", "rcn_150134", "rcn_150132", "rcn_150140", "rcn_150138", "rcn_14746", "rcn_14745", "rcn_14748", "rcn_14747", "rcn_14750", "rcn_150213", "rcn_150214", "rcn_198959", "rcn_14760rb", "rcn_14761"],
}

pca = ipa.pca(
data=data,
imap=imap,
minmap=0,
mincov=0.95,
impute_method="sample",
)

pca.run(nreplicates=100, seed=2146)
df = pd.DataFrame(pca.pcaxes[0], index=pca.names)
df.to_csv("pca_reconditus95.csv")
dfv = pd.DataFrame(pca.variances)
dfv.to_csv("var_reconditus95.csv")
dfm = pd.DataFrame(pca.missing)
dfm.to_csv("miss_reconditus95.csv")
'

sed -i "1s/,/,PC/g" pca_reconditus95.csv
sed -i "1s/^/ID/g" pca_reconditus95.csv
sed -i "1s/,/,rep/g" var_reconditus95.csv
sed -i "1s/^/PC/g" var_reconditus95.csv
sed -i "1s/^/ID/g" miss_reconditus95.csv



# grahami + garmani 
python -c'

import ipyrad.analysis as ipa
import pandas as pd
import toyplot

data = "/n/holyscratch01/losos_lab/Users/imaayan/all_ddrad/ipyrad/gg/gg_G95_outfiles/gg_G95.snps.hdf5"

imap={
"gg":["gmn_14593", "gmn_14588", "gmn_14594", "gmn_14582", "gmn_150352", "gmn_150354", "gmn_14638", "gmn_14643", "gmn_14641", "gmn_14646", "gmn_150392", "gmn_150380", "gmn_150711", "gmn_150720", "gmn_150623", "gmn_150054", "gmn_150089", "gmn_150087", "gmn_150400", "gmn_149886", "gmn_14606", "gmn_14597", "gmn_14599", "gmn_14609", "gmn_14604", "gmn_150451", "gmn_150449", "gmn_149937", "gmn_149858", "gmn_149875", "gmn_149876", "gmn_14711", "gmn_14718ra", "gmn_150678", "gmn_150684", "gmn_149952", "gmn_150655", "gmn_150626", "gmn_150761", "gmn_150759", "gmn_149905", "gmn_198954", "gmn_150294", "gmn_149923", "gmn_14686", "gmn_14696", "gmn_14694", "gmn_14693", "gmn_150605", "gmn_150036", "gmn_150035", "gmn_150365", "gmn_150795", "gmn_150794", "gmn_14547", "gmn_14551", "gmn_14763", "gmn_14757", "gmn_14755", "gmn_150536", "gmn_150555", "gmn_14568", "gmn_14563", "gmn_14569", "gmn_149990", "gmn_14618", "gmn_14623", "gmn_150315", "gmn_150311", "gmn_150673", "gmn_149951", "gmn_150112", "gmn_150114", "gmn_149962", "gmn_149998", "gmn_149991", "gmn_150585", "gmn_150578", "gmn_14673", "gmn_14672", "gmn_150735", "gmn_150752", "gmn_149843", "gmn_150725", "gmn_150697", "gmn_149846", "gmn_150811", "gmn_150807", "gmn_150708", "gmn_149881", "gmn_150472", "gmn_150468", "gmn_150595", "gmn_150599", "gmn_14832", "gmn_150129", "gmn_150126", "gmn_14819", "gmn_14816rb", "gmn_14818", "gmn_150120", "gmn_150166", "grh_14581", "grh_14591", "grh_14578", "grh_149864", "grh_149860", "grh_149958", "grh_14634", "grh_14655", "grh_14656", "grh_150384", "grh_150390", "grh_150714", "grh_150715", "grh_150613", "grh_150608", "grh_150615", "grh_198955", "grh_150057", "grh_150069", "grh_150060", "grh_150059", "grh_150077", "grh_150079", "grh_150086", "grh_150083", "grh_150799", "grh_150798", "grh_149894", "grh_149892", "grh_14614", "grh_14605", "grh_14602", "grh_150458", "grh_2674", "grh_2641", "grh_14784", "grh_14788", "grh_198958", "grh_150524", "grh_150521", "grh_149849", "grh_2581", "grh_2659", "grh_2658", "grh_2637", "grh_2663", "grh_2639", "grh_14710", "grh_14723", "grh_14722", "grh_14706", "grh_14715", "grh_149950", "grh_149943", "grh_150625", "grh_150650", "grh_3707", "grh_3683", "grh_3882ra", "grh_150765", "grh_150188", "grh_150192", "grh_150416", "grh_149909", "grh_150170", "grh_150169", "grh_150282", "grh_150290", "grh_150814", "grh_150482", "grh_150174", "grh_150179", "grh_150606", "grh_150037", "grh_2613", "grh_150260", "grh_150245", "grh_150154", "grh_150141", "grh_150783", "grh_14552", "grh_14727", "grh_150264", "grh_150261", "grh_14740", "grh_14741", "grh_14742", "grh_14744", "grh_2606", "grh_150441", "grh_150534", "grh_150537", "grh_14560", "grh_14561", "grh_14570", "grh_14559", "grh_14557", "grh_149818", "grh_149825", "grh_150496", "grh_150497", "grh_150377", "grh_14677", "grh_14658", "grh_150322", "grh_150301", "grh_14680", "grh_14679", "grh_150669", "grh_150668", "grh_2760", "grh_2747", "grh_14800", "grh_14794", "grh_150008", "grh_14808", "grh_14797", "grh_14796", "grh_150116", "grh_150109", "grh_150104", "grh_150101", "grh_150205", "grh_150203", "grh_149967", "grh_150098", "grh_150587", "grh_150583", "grh_150581", "grh_150586", "grh_198962", "grh_149868", "grh_149866", "grh_149869", "grh_150200", "grh_150193", "grh_150223", "grh_150229", "grh_150234", "grh_2622", "grh_2577", "grh_2624", "grh_2623", "grh_150746", "grh_150758", "grh_150753", "grh_150699", "grh_150694", "grh_150343", "grh_150340", "grh_150810", "grh_150808", "grh_14826", "grh_14824", "grh_14825", "grh_150427", "grh_149899", "grh_149896", "grh_149898", "grh_150687", "grh_150596", "grh_150591", "grh_150604", "grh_150168", "grh_150162", "grh_150163", "grh_14537", "grh_14536", "grh_14539", "grh_150575", "grh_150564"],
}

pca = ipa.pca(
data=data,
imap=imap,
minmap=0,
mincov=0.95,
impute_method="sample",
)

pca.run(nreplicates=100, seed=2147)
df = pd.DataFrame(pca.pcaxes[0], index=pca.names)
df.to_csv("pca_gg95.csv")
dfv = pd.DataFrame(pca.variances)
dfv.to_csv("var_gg95.csv")
dfm = pd.DataFrame(pca.missing)
dfm.to_csv("miss_gg95.csv")
'

sed -i "1s/,/,PC/g" pca_gg95.csv
sed -i "1s/^/ID/g" pca_gg95.csv
sed -i "1s/,/,rep/g" var_gg95.csv
sed -i "1s/^/PC/g" var_gg95.csv
sed -i "1s/^/ID/g" miss_gg95.csv


# lineatopus + reconditus
python -c'

import ipyrad.analysis as ipa
import pandas as pd
import toyplot

data = "/n/holyscratch01/losos_lab/Users/imaayan/all_ddrad/ipyrad/lr/lr_G95_outfiles/lr_G95.snps.hdf5"

imap={
"lr":["ltp_14589", "ltp_14592", "ltp_14576", "ltp_149863", "ltp_149859", "ltp_150331", "ltp_150330", "ltp_150329", "ltp_150353", "ltp_150361", "ltp_14657", "ltp_14654", "ltp_150382", "ltp_150387", "ltp_150381", "ltp_150716", "ltp_150712", "ltp_150619", "ltp_150617", "ltp_150610", "ltp_150064", "ltp_150068", "ltp_150070", "ltp_150072", "ltp_150393", "ltp_149888", "ltp_149893", "ltp_14603", "ltp_14600", "ltp_14608", "ltp_14607", "ltp_14782", "ltp_14786", "ltp_14787", "ltp_150527", "ltp_150529", "ltp_150530", "ltp_150514", "ltp_149850", "ltp_149851", "ltp_149856", "ltp_149857", "ltp_149874", "ltp_149878", "ltp_149872", "ltp_14709", "ltp_14705", "ltp_149947", "ltp_149946", "ltp_150682", "ltp_150628", "ltp_150643", "ltp_150640", "ltp_3829", "ltp_3469", "ltp_3494", "ltp_150767", "ltp_198964", "ltp_150187", "ltp_150189", "ltp_150407", "ltp_150413", "ltp_149920", "ltp_149914", "ltp_149907", "ltp_149977", "ltp_150284", "ltp_150289", "ltp_150479", "ltp_150487", "ltp_149929", "ltp_149930", "ltp_149925", "ltp_149931", "ltp_14692", "ltp_14688", "ltp_14684", "ltp_14683", "ltp_14682rb", "ltp_150175", "ltp_150172", "ltp_150032", "ltp_150030", "ltp_150374", "ltp_150248", "ltp_150240", "ltp_150258", "ltp_150152", "ltp_150150", "ltp_150780", "ltp_150797", "ltp_14548", "ltp_14543", "ltp_14544", "ltp_14550", "ltp_14540", "ltp_14731", "ltp_14733", "ltp_150263", "ltp_150265", "ltp_14736", "ltp_14739", "ltp_150436", "ltp_150431", "ltp_150539", "ltp_14572", "ltp_14558", "ltp_14566", "ltp_14562", "ltp_198966", "ltp_149822", "ltp_149820", "ltp_149831", "ltp_150501", "ltp_150499", "ltp_149989", "ltp_150375", "ltp_14629", "ltp_14628", "ltp_150324", "ltp_150327", "ltp_150305", "ltp_150325", "ltp_150303", "ltp_150304", "ltp_14676", "ltp_150666", "ltp_150670", "ltp_14804", "ltp_14801", "ltp_14805", "ltp_150018", "ltp_150019", "ltp_14810", "ltp_14809", "ltp_150110", "ltp_150115", "ltp_150111", "ltp_150100", "ltp_150105", "ltp_150106", "ltp_150212", "ltp_150206", "ltp_149961", "ltp_149993", "ltp_150006", "ltp_150094", "ltp_150090", "ltp_149971", "ltp_149969", "ltp_149972", "ltp_149975", "ltp_150049", "ltp_150052", "ltp_150042", "ltp_150198", "ltp_150196", "ltp_150239", "ltp_150221", "ltp_150222", "ltp_150238", "ltp_14670", "ltp_14671", "ltp_14663", "ltp_14668", "ltp_14665", "ltp_14669", "ltp_150425", "ltp_150426", "ltp_150417", "ltp_150428", "ltp_14773", "ltp_198957", "ltp_149982", "ltp_149985", "ltp_149988", "ltp_14765", "ltp_14772", "ltp_14768", "ltp_14769", "ltp_150736", "ltp_150739", "ltp_150751", "ltp_150754", "ltp_150728", "ltp_150724", "ltp_150704", "ltp_150696", "ltp_149842", "ltp_149840", "ltp_150339", "ltp_150342", "ltp_150806", "ltp_150812", "ltp_14822", "ltp_14823", "ltp_14821", "ltp_14820", "ltp_150707", "ltp_150709", "ltp_149882", "ltp_149880", "ltp_149883", "ltp_149901", "ltp_149904", "ltp_149902", "ltp_149903", "ltp_150469", "ltp_150465", "ltp_150462", "ltp_150685", "ltp_149953", "ltp_149955", "ltp_149954", "ltp_149956", "ltp_150602", "ltp_150600", "ltp_14836", "ltp_14834", "ltp_14833", "ltp_150117", "ltp_14815", "ltp_14817", "ltp_150122", "ltp_150118", "ltp_150160", "ltp_14530", "ltp_14534", "ltp_14522", "ltp_14523", "ltp_198961", "ltp_150570", "ltp_150566", "ltp_150560", "rcn_150133", "rcn_150136", "rcn_198956", "rcn_150139", "rcn_150134", "rcn_150132", "rcn_150140", "rcn_150138", "rcn_14746", "rcn_14745", "rcn_14748", "rcn_14747", "rcn_14750", "rcn_150213", "rcn_150214", "rcn_198959", "rcn_14760rb", "rcn_14761"],
}

pca = ipa.pca(
data=data,
imap=imap,
minmap=0,
mincov=0.95,
impute_method="sample",
)

pca.run(nreplicates=100, seed=2148)
df = pd.DataFrame(pca.pcaxes[0], index=pca.names)
df.to_csv("pca_lr95.csv")
dfv = pd.DataFrame(pca.variances)
dfv.to_csv("var_lr95.csv")
dfm = pd.DataFrame(pca.missing)
dfm.to_csv("miss_lr95.csv")
'

sed -i "1s/,/,PC/g" pca_lr95.csv
sed -i "1s/^/ID/g" pca_lr95.csv
sed -i "1s/,/,rep/g" var_lr95.csv
sed -i "1s/^/PC/g" var_lr95.csv
sed -i "1s/^/ID/g" miss_lr95.csv
