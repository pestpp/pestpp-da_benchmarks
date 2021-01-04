import os
import sys
import shutil
import platform
import numpy as np
import pandas as pd
import platform
import pyemu

bin_path = os.path.join("test_bin")
if "linux" in platform.platform().lower():
    bin_path = os.path.join(bin_path,"linux")
elif "darwin" in platform.platform().lower():
    bin_path = os.path.join(bin_path,"mac")
else:
    bin_path = os.path.join(bin_path,"win")

bin_path = os.path.abspath("test_bin")
os.environ["PATH"] += os.pathsep + bin_path


bin_path = os.path.join("..","..","..","bin")
exe = ""
if "windows" in platform.platform().lower():
    exe = ".exe"
exe_path = os.path.join(bin_path, "pestpp-ies" + exe)


noptmax = 4
num_reals = 20
port = 4021


def da_prep_4_mf6_freyberg_seq():
    t_d = os.path.join("mf6_freyberg","template_seq")
    if os.path.exists(t_d):
        shutil.rmtree(t_d)
    shutil.copytree(os.path.join("mf6_freyberg","template"),t_d)
    for f in os.listdir(t_d):
        for tag in ["ies","opt","glm"]:
            if tag in f.lower():
                os.remove(os.path.join(t_d,f))

    #first modify the tdis
    with open(os.path.join(t_d,"freyberg6.tdis"),'w') as f:
        f.write("BEGIN Options\n  TIME_UNITS  days\nEND Options\n")
        f.write("BEGIN Dimensions\n  NPER  1\nEND Dimensions\n")
        f.write("BEGIN PERIODDATA\n31.00000000  1       1.00000000\nEND PERIODDATA\n")
    #make sure it runs
    pyemu.os_utils.run("mf6",cwd=t_d)

    # write a tdis template file - could possibly keep all 25 stress periods to
    # simulate a 2-year-ahead forecast...

    with open(os.path.join(t_d,"freyberg6.tdis.tpl"),'w') as f:
        f.write("ptf  ~\n")
        f.write("BEGIN Options\n  TIME_UNITS  days\nEND Options\n")
        f.write("BEGIN Dimensions\n  NPER  1\nEND Dimensions\n")
        f.write("BEGIN PERIODDATA\n~  perlen  ~  1       1.00000000\nEND PERIODDATA\n")
    new_tpl,new_in = [os.path.join(t_d,"freyberg6.tdis.tpl")],[os.path.join(t_d,"freyberg6.tdis")]
    new_tpl_cycle = [-1]

    # mod the sto to make sp 1 transient - or should this be a pre-processor so that cycle 0 is 
    # ss and the rest are transient?

    # split out the head, sfr and list instruction files into multiple instruction file
    #eventually, want to move to da par and obs cycle tables for heads and gage_1 obs
    lines = open(os.path.join(t_d,"heads.csv.ins"),'r').readlines()[2:]
    new_ins,new_out,new_ins_cycle = [],[],[]
    #print(lines)
    for icycle, line in enumerate(lines):
        ins_name = os.path.join(t_d,"heads_{0}.csv.ins".format(icycle))
        with open(ins_name,'w') as f:
            f.write("pif ~\nl1\n")
            f.write(line)
        new_ins.append(ins_name)
        new_out.append(os.path.join(t_d,"heads.csv"))
        new_ins_cycle.append(icycle)
    remove_ins = ["heads.csv.ins"]

    lines = open(os.path.join(t_d,"sfr.csv.ins"),'r').readlines()[2:]
    #print(lines)
    for icycle, line in enumerate(lines):
        ins_name = os.path.join(t_d,"sfr_{0}.csv.ins".format(icycle))
        with open(ins_name,'w') as f:
            f.write("pif ~\nl1\n")
            f.write(line)
        new_ins.append(ins_name)
        new_out.append(os.path.join(t_d,"sfr.csv"))
        new_ins_cycle.append(icycle)
    remove_ins.append("sfr.csv.ins")

    lines = open(os.path.join(t_d,"freyberg6.lst.ins"),'r').readlines()[1:]
    icycle = 0
    tag_line = lines[0]
    for s in range(0,len(lines),13):
        ins_name = os.path.join(t_d,"freyberg6_{0}.lst.ins".format(icycle))
        with open(os.path.join(t_d,"freyberg6_{0}.lst.ins".format(icycle)),'w') as f:
            f.write("pif ~\n")
            f.write(tag_line)
            for line in lines[s+1:s+13]:
                f.write(line)
        new_ins.append(ins_name)
        new_out.append(os.path.join(t_d,"freyberg6.lst"))
        new_ins_cycle.append(icycle)
        icycle += 1
    remove_ins.append("freyberg6.lst.ins")

    # modify the ic file
    k = 0
    with open(os.path.join(t_d,"freyberg6.ic"),'r') as f:
        while True:
            line = f.readline()
            if line == "":
                break
            if line.lower().strip().startswith("internal"):
                arr_lines = []
                while True:
                    line = f.readline()
                    if line == "":
                        raise Exception
                    if line.lower().strip().startswith("end"):
                        break
                    if line.lower().strip().startswith("internal"):
                        arr = np.array(arr_lines,dtype=np.float)
                        np.savetxt(os.path.join(t_d,"heads_{0}.dat_in".format(k)),arr,fmt="%15.6E")
                        k += 1
                        arr_lines = []
                    else:
                        arr_lines.append(line.strip().split())

        arr = np.array(arr_lines, dtype=np.float)
        np.savetxt(os.path.join(t_d, "heads_{0}.dat_in".format(k)), arr, fmt="%15.6E")
    with open(os.path.join(t_d,"freyberg6.ic"),'w') as f:
        f.write("begin griddata\nstrt layered\n")
        for k in range(3):
            f.write("open/close 'heads_{0}.dat_in' FACTOR 1.0\n".format(k))
        f.write("end griddata\n")


    # write a python script to extract final heads and save to files
    with open(os.path.join(t_d,"forward_run.py"),'w') as f:
        f.write("import numpy as np\nimport flopy\nimport pyemu\n")
        f.write("pyemu.os_utils.run('mf6')\n")
        f.write("hds = flopy.utils.HeadFile('freyberg6_freyberg.hds')\n")
        f.write("arr = hds.get_data()\n")
        f.write("for k,a in enumerate(arr):\n")
        f.write("    np.savetxt('heads_'+str(k)+'.dat',a,fmt='%15.6E')\n")

    # dont run it so we can harvest the ic values in the arrays for setting the parval1 values
    pyemu.os_utils.run("python forward_run.py",cwd=t_d)

    # now write ins and tpl file for these
    ic_parvals = {}
    for k in range(3):
        fname = os.path.join(t_d,"heads_{0}.dat".format(k))
        assert os.path.exists(fname),fname
        arr = np.loadtxt(fname)
        fname_ins = fname + "_out.ins"
        fname_tpl = fname + "_in.tpl"
        in_arr = np.loadtxt(fname_tpl.replace(".tpl",""))
        ft = open(fname_tpl,'w')
        with open(fname_ins,'w') as f:
            f.write("pif ~\n")
            ft.write("ptf ~\n")
            for i in range(arr.shape[0]):
                f.write("l1 ")
                for j in range(arr.shape[1]):
                    if np.abs(arr[i,j]) > 100 or np.abs(in_arr[i,j]) > 100:
                        f.write(" !dum! ")
                        ft.write(" 40 ")
                    else:
                        oname = "head_{0:02d}_{1:03d}_{2:03d}".format(k,i,j)
                        f.write(" !{0}! ".format(oname))
                        ft.write(" ~  {0} ~ ".format(oname))
                        ic_parvals[oname] = in_arr[i,j]
                f.write("\n")
                ft.write("\n")
        ft.close()
        new_tpl.append(fname_tpl)
        new_in.append(fname_tpl.replace(".tpl",""))
        new_tpl_cycle.append(-1)
        new_ins.append(fname_ins)
        new_out.append(os.path.join(t_d,"heads_{0}.dat".format(k)))
        new_ins_cycle.append(-1)

        i = pyemu.pst_utils.InstructionFile(fname_ins)
        df = i.read_output_file(fname)
        #print(df)

    # split out the wel and rch tpl files into cycle files
    lines = []
    with open(os.path.join(t_d,"freyberg6.wel.tpl"),'r') as f:
        for i in range(19):
            lines.append(f.readline())
    print(lines)

    for icycle in range(25):
        tpl_file = os.path.join(t_d,"freyberg6.wel_{0}.tpl".format(icycle))
        with open(tpl_file,'w') as f:
            for line in lines:
                new_line = line.replace("_0","_{0}".format(icycle))
                f.write(new_line)

        new_tpl.append(tpl_file)
        new_in.append(os.path.join(t_d,"freyberg6.wel"))
        new_tpl_cycle.append(icycle)
    remove_tpl = ["freyberg6.wel.tpl"]

    lines = []
    with open(os.path.join(t_d, "freyberg6.rch.tpl"), 'r') as f:
        for i in range(11):
            lines.append(f.readline())
    print(lines)

    for icycle in range(25):
        tpl_file = os.path.join(t_d,"freyberg6.rch_{0}.tpl".format(icycle))
        with open(tpl_file,'w') as f:
            for line in lines:
                new_line = line.replace("_0","_{0}".format(icycle))
                f.write(new_line)

        new_tpl.append(tpl_file)
        new_in.append(os.path.join(t_d,"freyberg6.rch"))
        new_tpl_cycle.append(icycle)
    remove_tpl.append('freyberg6.rch.tpl')

    # now for the fun part: modify the pst
    pst = pyemu.Pst(os.path.join(t_d,"freyberg6_run.pst"))

    print(pst.npar_adj,pst.nnz_obs)

    # swap out obs info
    dropped_dfs = []
    for ins in remove_ins:
        dropped_dfs.append(pst.drop_observations(os.path.join(t_d, ins), '.'))
    for insf, outf, cy in zip(new_ins, new_out, new_ins_cycle):
        df = pst.add_observations(insf,outf, pst_path=".")
        pst.observation_data.loc[df.obsnme, "cycle"] = cy
        pst.model_output_data.loc[os.path.join(".",os.path.split(insf)[1]),"cycle"] = cy
    pst.observation_data.loc[:,"weight"] = 0.0
    for df in dropped_dfs:
        for c in ["obsval","weight"]:
            pst.observation_data.loc[df.obsnme, c] = df.loc[:, c]

    # swap out par info
    dropped_dfs = []
    for tpl in remove_tpl:
        dropped_dfs.append(pst.drop_parameters(os.path.join(t_d,tpl),'.'))
    pst.parameter_data.loc[:,"cycle"] = -1
    pst.model_input_data.loc[:,"cycle"] = -1
    for tplf, inf, cy in zip(new_tpl,new_in,new_tpl_cycle):
        df = pst.add_parameters(tplf,inf,pst_path=".")
        pst.parameter_data.loc[df.parnme,"cycle"] = cy
        pst.model_input_data.loc[os.path.join(".",os.path.split(tplf)[1]),"cycle"] = cy
    for df in dropped_dfs:
        for c in ["parval1","parubnd","parlbnd","pargp"]:
            pst.parameter_data.loc[df.parnme,c] = df.loc[:,c]

    #set the state parameter info
    for p,v in ic_parvals.items():
        pst.parameter_data.loc[p,"parval1"] = v
        pst.parameter_data.loc[p, "parlbnd"] = v * 0.5
        pst.parameter_data.loc[p, "parubnd"] = v * 1.5
        pst.parameter_data.loc[p,"pargp"] = "head_state"
        pst.parameter_data.loc[p,"partrans"] = "none"

    pst.control_data.noptmax = 2
    pst.model_command = "python forward_run.py"
    pst.pestpp_options.pop("ies_par_en")
    pst.parameter_data.loc["perlen","partrans"] = "fixed"
    pst.pestpp_options["ies_num_reals"] = 5
    pst.pestpp_options["da_num_reals"] = 5
    pst.write(os.path.join(t_d,"freyberg6_run_da1.pst"),version=2)
    return pst

def da_mf6_freyberg_test_1():
    test_d = "mf6_freyberg"
    t_d = os.path.join(test_d,"template_seq")
    
    pst = da_prep_4_mf6_freyberg_seq()
    
    print(exe_path.replace("ies","da"))
    pst.control_data.noptmax = 0
    pst.pestpp_options["ies_verbose_level"] = 4
    pst.pestpp_options["ies_no_noise"] = True
    pst.write(os.path.join(t_d, "freyberg6_run_da1.pst"), version=2)
    pyemu.os_utils.run("{0} freyberg6_run_da1.pst".format(exe_path.replace("ies","da")),cwd=t_d)

    pst.pestpp_options["ies_num_reals"] = 15
    pst.control_data.noptmax = 2
    pst.write(os.path.join(t_d, "freyberg6_run_da1.pst"), version=2)
    pyemu.os_utils.start_workers(t_d,exe_path.replace("ies","da"),"freyberg6_run_da1.pst",
                                 num_workers=15,worker_root=test_d,port=port,
                                 master_dir=os.path.join(test_d,"master_da_1"),verbose=True)
    
def da_prep_4_mf6_freyberg_seq_tbl():
    test_d = "mf6_freyberg"
    t_d = os.path.join(test_d,"template_seq")
    if not os.path.exists(t_d):
        da_prep_4_mf6_freyberg_seq()
    pst = pyemu.Pst(os.path.join(t_d,"freyberg6_run_da1.pst"))
    pdf = pd.DataFrame({"perlen":31},index=np.arange(25))
    pdf.T.to_csv(os.path.join(t_d,"par_cycle_tbl.csv"))
    pst.pestpp_options["da_parameter_cycle_table"] = "par_cycle_tbl.csv"

    # mod sfr_out
    sfr_ins_file = os.path.join(t_d, "sfr.csv.ins")
    with open(sfr_ins_file, 'w') as f:
        f.write("pif ~\n")
        f.write("l1\n")
        f.write("l1 ~,~ !headwater!  ~,~ !tailwater!  ~,~ !gage_1!\n")

    # and lst budget
    lines = open(os.path.join(t_d, "freyberg6_0.lst.ins"), 'r').readlines()
    lst_ins_file = os.path.join(t_d, "freyberg6.lst.ins")
    with open(lst_ins_file, 'w') as f:
        for line in lines:
            f.write(line.replace("_20151231",""))

    obs = pst.observation_data.copy()
    tr_obs = obs.loc[obs.obsnme.str.startswith("trgw"),:].copy()
    tr_obs.loc[tr_obs.obsnme,"datetime"] = pd.to_datetime(tr_obs.obsnme.apply(lambda x: x.split('_')[-1]))
    tr_obs.loc[tr_obs.obsnme,"k"] = tr_obs.obsnme.apply(lambda x: np.int(x.split('_')[1]))
    tr_obs.loc[tr_obs.obsnme, "i"] = tr_obs.obsnme.apply(lambda x: np.int(x.split('_')[2]))
    tr_obs.loc[tr_obs.obsnme, "j"] = tr_obs.obsnme.apply(lambda x: np.int(x.split('_')[3]))
    tr_obs.loc[tr_obs.obsnme,"obgnme"] = tr_obs.obsnme.apply(lambda x: "_".join(x.split("_")[:-1]))

    head_obs = obs.loc[obs.obsnme.str.startswith("head_"),:].copy()
    head_obs.loc[head_obs.obsnme, "k"] = head_obs.obsnme.apply(lambda x: np.int(x.split('_')[1]))
    head_obs.loc[head_obs.obsnme, "i"] = head_obs.obsnme.apply(lambda x: np.int(x.split('_')[2]))
    head_obs.loc[head_obs.obsnme, "j"] = head_obs.obsnme.apply(lambda x: np.int(x.split('_')[3]))

    #print(pst.nobs)
    for ins_file in pst.model_output_data.pest_file.copy():
        if "heads_" in ins_file and ins_file.endswith("dat_out.ins"):
            continue
        pst.drop_observations(os.path.join(t_d,ins_file),pst_path=".")
    for ins_file in [sfr_ins_file,lst_ins_file]:
        pst.add_observations(ins_file,pst_path=".")

    # work out which heads are obs site
    obs_heads = {}
    odf_names = []
    pst.observation_data.loc[:,"org_obgnme"] = np.NaN
    pst.observation_data.loc[:, "weight"] = 0.0
    for og in tr_obs.obgnme.unique():
        site_obs = tr_obs.loc[tr_obs.obgnme==og,:]
        site_obs.sort_values(by="datetime",inplace=True)
        head_name = "head_{0:02d}_{1:03d}_{2:03d}".format(site_obs.k[0],site_obs.i[0],site_obs.j[0])
        for i,oname in enumerate(site_obs.obsnme):
            obs_heads[oname] = (head_name,i)
        # assign max weight in the control file since some have zero weight and
        # we are covering weights in the weight table
        #print(og,site_obs.weight.max())
        pst.observation_data.loc[head_name,"weight"] = site_obs.weight.max()
        pst.observation_data.loc[head_name,"org_obgnme"] = og
        odf_names.append(head_name)

    odf_names.append("gage_1")

    odf = pd.DataFrame(columns=odf_names,index=np.arange(25))
    wdf = pd.DataFrame(columns=odf_names,index=np.arange(25))
    for tr_name,(head_name,icycle) in obs_heads.items():
        odf.loc[icycle,head_name] = obs.loc[tr_name,"obsval"]
        wdf.loc[icycle, head_name] = obs.loc[tr_name, "weight"]

    g_obs = obs.loc[obs.obsnme.str.startswith("gage_1"),:].copy()
    #give these obs the max weight since some have zero weight
    pst.observation_data.loc["gage_1", "weight"] = g_obs.weight.max()
    g_obs.sort_index(inplace=True)
    for i,name in enumerate(g_obs.obsnme):
        odf.loc[i,"gage_1"] = g_obs.loc[name,"obsval"]
        wdf.loc[i, "gage_1"] = g_obs.loc[name, "weight"]

    # now drop any entries that have zero weight across all cycles
    print(pst.nnz_obs_names)
    odf = odf.loc[:,pst.nnz_obs_names]
    wdf = wdf.loc[:, pst.nnz_obs_names]

    odf.T.to_csv(os.path.join(t_d,"obs_cycle_tbl.csv"))
    pst.pestpp_options["da_observation_cycle_table"] = "obs_cycle_tbl.csv"
    wdf.T.to_csv(os.path.join(t_d, "weight_cycle_tbl.csv"))
    pst.pestpp_options["da_weight_cycle_table"] = "weight_cycle_tbl.csv"

    pst.observation_data.loc[:,"cycle"] = -1
    pst.model_output_data.loc[:,"cycle"] = -1

    pst.write(os.path.join(t_d,"freyberg6_run_da2.pst"),version=2)
    return pst

def da_mf6_freyberg_test_2():
    pst = da_prep_4_mf6_freyberg_seq_tbl()
    test_d = "mf6_freyberg"
    t_d = os.path.join(test_d, "template_seq")
    print(exe_path.replace("ies", "da"))
    pst.control_data.noptmax = 0
    pst.pestpp_options["ies_verbose_level"] = 4
    pst.pestpp_options["ies_no_noise"] = True
    pst.write(os.path.join(t_d, "freyberg6_run_da2.pst"), version=2)
    pyemu.os_utils.run("{0} freyberg6_run_da2.pst".format(exe_path.replace("ies","da")),cwd=t_d)

    pst.control_data.noptmax = -1
    pst.pestpp_options["ies_num_reals"] = 15
    pst.write(os.path.join(t_d, "freyberg6_run_da2.pst"), version=2)
    pyemu.os_utils.start_workers(t_d, exe_path.replace("ies", "da"), "freyberg6_run_da2.pst",
                                num_workers=15, worker_root=test_d, port=port,
                                master_dir=os.path.join(test_d, "master_da_2"), verbose=True)


def invest():
    t_d = os.path.join("mf6_freyberg","template_seq")
    pst1 = pyemu.Pst(os.path.join(t_d,"freyberg6_run_da1.pst"))
    pst2 = pyemu.Pst(os.path.join(t_d,"freyberg6_run_da2.pst"))
    print(pst1.observation_data.loc[pst1.nnz_obs_names,["obsval","weight"]])
    odf = pd.read_csv(os.path.join(t_d,pst2.pestpp_options["da_observation_cycle_table"]))
    wdf = pd.read_csv(os.path.join(t_d, pst2.pestpp_options["da_weight_cycle_table"]))
    print(odf)
    print(wdf)

def da_mf6_freyberg_smoother_test():
    model_d = "mf6_freyberg"
    local=True
    if "linux" in platform.platform().lower() and "10par" in model_d:
        #print("travis_prep")
        #prep_for_travis(model_d)
        local=False
    
    t_d = os.path.join(model_d,"template")

    m_d = os.path.join(model_d,"master_da_smoother")
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    pst = pyemu.Pst(os.path.join(t_d,"freyberg6_run_ies.pst"))
    pst.control_data.noptmax = 0
    pst.pestpp_options["ies_autoadaloc"] = False
    pst.pestpp_options.pop("ies_localizer",None)
    pst.pestpp_options["da_num_reals"] = 15
    pst.write(os.path.join(t_d,"freyberg6_run_da.pst"))
    pyemu.os_utils.run("{0} freyberg6_run_da.pst".format(exe_path.replace("-ies","-da")),cwd=t_d)
    
    pst.control_data.noptmax = 2
    
    pst.write(os.path.join(t_d,"freyberg6_run_da.pst"))
    pyemu.os_utils.start_workers(t_d, exe_path.replace("-ies","-da"), "freyberg6_run_da.pst", num_workers=15,
                                master_dir=m_d,worker_root=model_d,port=port)

    
    oe_file = os.path.join(m_d,"freyberg6_run_da.global.0.oe.csv")
    assert os.path.exists(oe_file),oe_file
    pe_file = os.path.join(m_d,"freyberg6_run_da.global.prior.pe.csv")
    assert os.path.exists(pe_file),pe_file
    phi_file = os.path.join(m_d,"freyberg6_run_da.global.phi.actual.csv")
    assert os.path.exists(phi_file),phi_file


def build_lorenz_pst():
    t_d = os.path.join("lorenz96", "template_seq")
    if os.path.exists(t_d):
        shutil.rmtree(t_d)
    shutil.copytree(os.path.join("lorenz96", "template"), t_d)
    os.chdir(os.path.join(t_d))

    # run forward check
    pyemu.os_utils.run("python rk4_lorenz96.py")

    out_csv = "l96_out.csv"
    odf = pyemu.pst_utils.csv_to_ins_file(os.path.join(out_csv))# only_cols=cols)
    ins_file = out_csv + ".ins"

    in_csv = "l96_in.csv"
    F = pd.read_csv(os.path.join(in_csv))['F'][0]
    tpl_file = in_csv + ".tpl"

    pst = pyemu.helpers.pst_from_io_files(tpl_file, in_csv, ins_file, out_csv)

    obs = pst.observation_data
    obs.loc[odf.index, "obsval"] = odf.obsval

    par = pst.parameter_data
    par.loc["forcing", "parval1"] = F
    par.loc["forcing", "parlbnd"] = 0.0
    par.loc["forcing", "parubnd"] = 9.0

    pst.control_data.noptmax = 1
    pst.model_command = "python forward_run.py"

    with open(os.path.join("forward_run.py"), 'w') as f:
        f.write("import rk4_lorenz96\n")
        f.write("rk4_lorenz96.L96()\n")

    pst.write(pst.filename)

    #do pestpp test run here

    os.chdir(os.path.join("..", ".."))


if __name__ == "__main__":
    
    
    shutil.copy2(os.path.join("..","exe","windows","x64","Debug","pestpp-da.exe"),os.path.join("..","bin","pestpp-da.exe"))
    da_mf6_freyberg_test_1()
