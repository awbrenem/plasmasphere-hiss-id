#! /usr/bin/env python

#Loads the file containing all the entry and exit times for the plasmasphere and
#passes these to hiss_build_database_driver.py for each day from d0 to d1


import datetime
import time
import subprocess


#Select input parameters
p = 'a'
d0 = datetime.datetime(2014, 1, 3, 0, 0, 0)
d1 = datetime.datetime(2014, 1, 6, 0, 0, 0)


#Read file with all entry and exit times in plasmasphere
fn = 'plasmasphere_rbsp'+p+'_database.txt'
f = open(fn, 'r')
print f
vals = f.read()
vals = vals.split("\n")
type(vals)
len(vals)


strstart = []
strend = []
fullstart = []
fullend = []
for x in range(0,len(vals)-1):
    vtmp = vals[x]
    strstart.append(vtmp[0:10])
    strend.append(vtmp[20:30])
    fullstart.append(vtmp[0:19])
    fullend.append(vtmp[20:39])
    f.close()



    ndays = str(d1 - d0)

    start = 0
    end = ndays.find('days')
    ndays = ndays[start:end-1]

    dnew = d0.isoformat()
    end = dnew.find('T')
    dnew = dnew[0:end]

    
    for x in range(0,int(ndays)):

        goostart = []
        gooend = []

        print "current date = ",dnew

        #find all instances of current date
        for i,j in enumerate(strstart):
            if j == dnew:
                print "\t",vals[i]
                goostart.append(fullstart[i])
                gooend.append(fullend[i])

                for i in range(0,len(goostart)-1):
                    print "\t\t",goostart[i]

                    print "----"
                    for i in range(0,len(gooend)-1):
                        print "\t\t",gooend[i]




    exit_code = subprocess.call(['/Applications/itt/idl71/bin/idl','-e','hiss_build_database_rbsp','-args',
                                 '%s'%dnew,'%s'%p,'%s'%goostart,'%s'%gooend])


    tnew = time.mktime(d0.timetuple()) + 86400*(x+1)

    dnew = d0.fromtimestamp(tnew)
    dnew = dnew.isoformat()

    goo = dnew.find('T')
    dnew = dnew[0:goo]
