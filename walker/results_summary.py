#!/usr/bin/python



## ip_prob_q1350sq292_m3.lp.flow.gz
#[((3 142727718) (3 153747870) 0 3 1308062 5590200) ((3 153747870) (3 160865436) 0 0 904916 3852464) ((3 160865436) (3 169827113) 0 3 968130 4105772) ((19 19319387) (19 20277334) 0 0 120786 535195) ((19 20277334) (19 20297843) 0 6 3490 17428) ((19 20297843) (19 23618913) 0 0 404676 1553278) ((19 23618913) (19 23705088) 0 3 11140 48104) ((19 23705088) (19 23730556) 0 6 3760 19438) ((19 23730556) (19 23735924) 0 3 812 3110) ((23 60040) (23 488813) 0 6 32506 187423) ((23 488813) (23 2485209) 0 0 196816 774180) ((23 2485209) (23 3695988) 0 6 150448 885304) ((23 3695988) (23 6393503) 0 0 304366 1095160) ((23 6393503) (23 8289236) 0 6 233070 1066650) ((23 8289236) (23 8365651) 0 0 9998 28026) ((23 18537397) (23 22023147) 0 0 430410 1592856) ((23 22023147) (23 22053794) 0 6 4196 22436) ((23 22053794) (23 31527675) 0 0 1107242 3999970) ((23 31527675) (23 32840998) 0 6 144856 682236) ((23 32840998) (23 33684998) 0 0 91538 254194) ((23 33684998) (23 34921309) 0 6 134706 601206) ((23 34921309) (23 43296155) 0 0 1013828 3674522) ((23 43296155) (23 44512249) 0 6 151644 777782) ((23 44512249) (23 46280140) 0 0 215508 736796) ((23 46280140) (23 46309332) 0 6 4076 21262) ((23 46309332) (23 46313490) 0 0 494 1636) ((23 46313490) (23 46314793) 0 6 154 768) ((23 46314793) (23 57428925) 0 0 1238400 4266240) ((23 57428925) (23 57452965) 0 6 2958 15198) ((23 57452965) (23 58579928) 0 0 147720 543750)]
#[((23 6393503) (23 22023147) 3 6 2 38) ((23 46309332) (23 46314793) 2 6 2 49) ((23 46280140) (23 46313490) 3 6 2 48) ((23 8289236) (23 33684998) 0 6 2 50) ((23 57428925) (23 57452965) 1 6 2 31) ((19 23730556) (19 23735924) 2 3 2 56) ((23 2485209) (23 32840998) 1 6 2 40) ((19 23618913) (19 23705088) 3 3 2 25) ((19 20277334) (19 20297843) 1 6 2 41) ((23 43296155) (23 44512249) 1 6 2 43) ((23 488813) (23 3695988) 2 6 2 29)]


import sys


if len(sys.argv)!=2:
	print "%s file" % sys.argv[0]
	sys.exit(1)


h=open(sys.argv[1])
lines=h.readlines()
if len(lines)!=3:
	print >> sys.stderr, "WARNING WEIRD NUMBER OF LINES!"


a=eval(lines[1].replace(' ',',').strip())
b=eval(lines[2].replace(' ',',').strip())

print len(a),len(b),float(len(a))/(len(b)+1)

h.close()