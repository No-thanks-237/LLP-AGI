
import sys, os
RUN = r"E:\AAAwork\LLP-AGI\runs\bc5_lhcb_downstream_v3"
sys.path.insert(0, os.path.join(RUN, "work", "analysis"))
import bc5_lib as L
print("%7s %7s %7s %7s %7s %7s %7s %7s" % ("m","pipi","KK","multi","gg","ss","cc","tau"))
for m in [1.0,1.5,1.8,1.99,2.0,2.2,2.5,3.0]:
    br,_ = L.channel_branching(m)
    print("%7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f" % (
        m, br['pipi'], br['KK'], br['multi'], br['gg'], br['ss'], br['cc'], br['tau']))
