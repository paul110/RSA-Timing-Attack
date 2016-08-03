import sys, subprocess, hashlib, shlex, math, binascii, random, copy
import threading
from multiprocessing import Process


global ro
global ro_s
global omega

global sampleSize
sampleSize = 3500

global base
base = 2**64

global maxKeySize
maxKeySize = 256

global treshold
treshold = 15

global lookAhead
lookAhead = 3

global interactions
interactions = 0



def getLimb(x, i):
    global base
    result = (x >> 64*i) & (base-1)
    return result

def calc_ro(n):
    ro = 1
    while ro < n:
        ro = ro << 64
    return  ro

def interact(  c ) :
  global interactions
  interactions += 1
  target_in.write("{0:X}".format(c) + "\n")
  # target_in.write(n + "\n"  )
  # target_in.write(d  + "\n"  )

  time = int( target_out.readline().strip())
  message = int(target_out.readline().strip(), 16)

  return (time, message)

def limbsNr(x):
    return int(math.ceil(math.log( x, base )))

def calc_ro_square(n):
    bits = 2 * limbsNr(n) * 64
    res = 1
    for i in xrange(0, bits):
        res = res * 2
        while(res > n):
            res = res - n
    return res

def getOmega(n):
    w1 = - getLimb(n, 0)
    w = 1
    for i in xrange(0, 64):
        w = (w*w * w1) % base
    return w

def MontMul(x, y, n, omega):
    extra = 0
    w = 64
    mods = limbsNr(n)
    r = 0
    for i in xrange( 0, mods):
        yil = y & (base-1)
        y = y / 2**64
        u = (omega * (yil * (x & (base-1))+ r&(base-1) )) & (base-1)
        temp = x * yil
        temp2 = n * u

        r = temp + r + temp2
        r = r >> w
    if(r > n):
        extra = 1
        r = r -n
    return r, extra

def computeSamples(n):
    ro_s = calc_ro_square(n)
    omega = getOmega(n)

    messages = []
    mess_origin = []
    timings = []
    results = []
    for j in range(0, sampleSize):
        m = random.getrandbits(n.bit_length())
        while m > n :
            m = random.getrandbits(n.bit_length())
        mess_origin.append(m)
        messages.append(MontMul(mess_origin[j], ro_s, n, omega)[0])

        res = interact(m)
        timings.append(res[0])
        results.append(res[1])
    return (mess_origin, messages, timings, results)

def updateMessages(messages ,d, n):
    currentSet = messages[:]
    for i in range(0, len(messages)):
        currentSet[i] = MontMul(currentSet[i], currentSet[i], n, omega)[0]

    for j in range(1, len(d)):
        for i in range(0, len(messages)):
            if(d[j] == "1"):
                currentSet[i] = MontMul(currentSet[i], messages[i], n, omega)[0]
            currentSet[i] = MontMul(currentSet[i], currentSet[i], n , omega)[0]
    return currentSet

def checkKey( message, result, d):
    d1 = int(d + "1", 2)
    d0 = int(d + "0", 2)

    # choose which prediction was right and output key
    if pow(message, d1, n) == result:
        # print "key : ", "{0:X}".format(d1)
        return (1, "1")
    else :
        if pow(message, d0, n)== result:
            # print "key : ", "{0:X}".format(d0)
            return (1, "0")
    # print "no key found", "{0:X}".format(int(d, 2))
    return (0, "0")

def getKeySize(n):
    omega = getOmega(n)
    ro = calc_ro(n)
    print "Getting key SIZE"

    message = MontMul(ro, 1, n, omega)[0]

    res = interact(MontMul(1, 1, n, omega)[0])
    print res
    print MontMul(1, 1, n, omega)[0]
    print "key size recovery finished"
    # return keySize

class sampleMessages (threading.Thread):
    def __init__(self, n, sampleSize):
        threading.Thread.__init__(self)
        self.n = n
        self.sampleSize = sampleSize
    def run(self):
        ro_s = calc_ro_square(self.n)
        omega = getOmega(self.n)
        mess_origin = []
        while len(mess_origin) < self.sampleSize:
            m = random.getrandbits(n.bit_length())
            while m > n :
                m = random.getrandbits(n.bit_length())
            mess_origin.append(m)
        self.mess_origin = mess_origin
    def join(self):
        return self.mess_origin

def checkAnormal(diff1, diff2):
    if (abs(diff1 - diff2) < treshold) | ((diff1 < 0) & (diff2 < 0)):
        return 1
    return 0

def attack(n, e):
    global sampleSize

    # set start conditions
    valid = 0
    keySize = 64
    stableKey = "1"

    while (keySize < maxKeySize) &  (not valid):
        # re-initialise everything
        mess_origin = []
        messages    = []
        currentSet  = []
        results     = []
        time        = []
        stablekeySet= 0
        # generate sample messages of n bits
        generateMessages = sampleMessages(n, sampleSize)
        generateMessages.run()
        mess_origin += generateMessages.join()
        currentSetTest = []

        # set key to stable key
        d = stableKey
        # compute times taken to decrypt each message
        for i in range(0, len(mess_origin)):
            # get times for each message
            res = interact(mess_origin[i])
            results.append(res[1])
            time.append( res[0] )
            # transform each message into montgomery form
            messages.append(MontMul(mess_origin[i], ro_s, n, omega)[0])

        print "Sample Size: ", sampleSize
        print "Interactions: ", interactions
        print "Start from key bit ", len(d)

        currentSet = updateMessages(messages, d, n)
        set1 = currentSet[:]
        set0 = currentSet[:]
        stableSet = currentSet[:]

        #reset errors
        error = 0

        # discover the next bits untill the last bit which needs to be guessed
        while ((not valid) & (len(d) <= keySize) & (error < 7)):
            warning = 0
            j = len(d)
            groupTime = [0, 0, 0, 0, 0]
            groupSize = [0, 0, 0, 0, 0]

            for i in range(0, len(messages)):
                # assume bit j is 0
                (set0[i], extra) = MontMul(currentSet[i], currentSet[i], n, omega)
                # if extra reduction
                if extra :
                    groupTime[3] += time[i]
                    groupSize[3] +=1
                # if no extra reduction
                else :
                    groupTime[4] += time[i]
                    groupSize[4] += 1

                # assume bit j is 1
                (temp, extra) = MontMul(currentSet[i], messages[i], n , omega)
                (set1[i], extra) = MontMul(temp, temp, n ,omega)
                # if extra reduction
                if extra :
                    groupTime[1] += time[i]
                    groupSize[1] +=1
                # if no extra reduction
                else :
                    groupTime[2] += time[i]
                    groupSize[2] +=1


            # compute averaget time for all 4 groups
            uF1 = float(groupTime[1])/groupSize[1]
            uF2 = float(groupTime[2])/groupSize[2]
            uF3 = float(groupTime[3])/groupSize[3]
            uF4 = float(groupTime[4])/groupSize[4]

            # compute differences between pari groups
            diff1 = uF1 - uF2
            diff2 = uF3 - uF4

            # condition for anormal behaviour
            if checkAnormal(diff1, diff2) :
                warning = 1
                error += 2
            else :
                if error > 0:
                    error -=1
                warning = 0

            if warning > 0 :
                # print "--------------------DEJA EROARE ",error, " ORI---------------"
                print "warning at bit", len(d)
                stablekeySet = 1
                # branch to 0 and 1
                saveSet0 = set0[:]
                saveSet1 = set1[:]
                turnsDelta = [0, 0]
                flags = [0, 0]
                currentSet0  = set0[:]
                currentSet1 = set1[:]
                d1 = d + "1"
                d0 = d + "0"
                (valid, bit) = checkKey(mess_origin[0], results[0], d0)
                if valid :
                    d = d0 + str(bit)
                    return d
                    # print "!!!!!!!!!!!!!!FoundKey!!!!!!!!!!!!!!!"
                else :
                    (valid, bit) = checkKey(mess_origin[0], results[0], d1)
                    if valid :
                        d = d1 + str(bit)
                        return d
                        # print "!!!!!!!!!!!!!!FoundKey!!!!!!!!!!!!!!!"
                if not valid :
                    # print "ver1", len(d0), ":",  d0[len(d0)-1], "----------------", int(diff1), int(diff2)

                    # check following rounds for bit 0
                    for y in range(0, lookAhead):
                        groupTime = [0, 0, 0, 0, 0]
                        groupSize = [0, 0, 0, 0, 0]

                        for k in range(0, len(currentSet)):

                            (set0[k], extra) = MontMul(currentSet0[k], currentSet0[k], n, omega)
                            if extra :
                                groupTime[3] += time[k]
                                groupSize[3] += 1
                            else :
                                groupTime[4] += time[k]
                                groupSize[4] += 1

                            (temp, extra) = MontMul(currentSet0[k], messages[k], n, omega)
                            (set1[k], extra) = MontMul(temp, temp, n, omega)
                            if extra :
                                groupTime[1] += time[k]
                                groupSize[1] +=1
                            # if no extra reduction
                            else :
                                groupTime[2] += time[k]
                                groupSize[2] +=1

                        uF1 = float(groupTime[1])/groupSize[1]
                        uF2 = float(groupTime[2])/groupSize[2]
                        uF3 = float(groupTime[3])/groupSize[3]
                        uF4 = float(groupTime[4])/groupSize[4]

                        diff1 = uF1 - uF2
                        diff2 = uF3 - uF4

                        bit = 0
                        if diff1 > diff2:
                            bit = 1
                        if checkAnormal(diff1, diff2):
                            flags[0] = 1

                        if bit :
                            currentSet0 = set1[:]
                        else:
                            currentSet0 = set0[:]
                        d0 += str(bit)
                        # print "ver1 : ", len(d0), ": ",  bit, "------", abs(int(diff1 - diff2)), "--------", int(diff1), int(diff2), "error:", error
                        (valid, bit) = checkKey(mess_origin[0], results[0], d0)
                        if valid :
                            d = d0 + str(bit)
                            return d
                            # print "!!!!!!!!!!!!!!FoundKey!!!!!!!!!!!!!!!"

                        if not checkAnormal(diff1, diff2) :
                            turnsDelta[0] += abs(int(diff1 - diff2))


                    # print "ver2", len(d1), ":",  d1[len(d1)-1], "---------------"
                    # check following rounds for bit 0
                    for y in range(0, lookAhead):
                        groupTime = [0, 0, 0, 0, 0]
                        groupSize = [0, 0, 0, 0, 0]
                        for k in range(0, len(currentSet1)):

                            (set0[k], extra) = MontMul(currentSet1[k], currentSet1[k], n, omega)
                            if extra :
                                groupTime[3] += time[k]
                                groupSize[3] += 1
                            else :
                                groupTime[4] += time[k]
                                groupSize[4] += 1

                            (temp, extra) = MontMul(currentSet1[k], messages[k], n, omega)
                            (set1[k], extra) = MontMul(temp, temp, n, omega)
                            if extra :
                                groupTime[1] += time[k]
                                groupSize[1] +=1
                            # if no extra reduction
                            else :
                                groupTime[2] += time[k]
                                groupSize[2] +=1

                        uF1 = float(groupTime[1])/groupSize[1]
                        uF2 = float(groupTime[2])/groupSize[2]
                        uF3 = float(groupTime[3])/groupSize[3]
                        uF4 = float(groupTime[4])/groupSize[4]

                        diff1 = uF1 - uF2
                        diff2 = uF3 - uF4

                        bit = 0
                        if diff1 > diff2:
                            bit = 1
                        if checkAnormal(diff1, diff2):
                            flags[1] =  1
                        if bit :
                            currentSet1 = set1[:]
                        else:
                            currentSet1 = set0[:]
                        d1 += str(bit)
                        # print "ver2:", len(d1), ": ",  bit, "------", abs(int(diff1 - diff2)), "--------", int(diff1), int(diff2), "error:", error
                        # if( diff1 == 0) :
                            # print "DIFF1 = 0", uF1, uF2
                        (valid, bit) = checkKey(mess_origin[0], results[0], d1)
                        if valid :
                            d = d1 + str(bit)
                            return d
                        if not ((diff1 < 0 ) & (diff2 < 0)):
                            turnsDelta[1] += abs(int(diff1 - diff2))

                    # decide which set was better
                    set0 = saveSet0[:]
                    set1 = saveSet1[:]
                    if(turnsDelta[0] < turnsDelta[1]) :
                        bit = 1
                        # print "CHOOSE 2"
                    else :
                        # print "CHOOSE 1"
                        bit = 0
            else :
                bit = 0
                #if diff for the groups with assumed bit = 1 is bigger than the diff of the groups with assumed bit = 0, then predict 1, else predict 0
                if diff1 > diff2:
                    bit = 1
                if(not stablekeySet):
                    stableKey = d

            if not valid :
                # depending on which bit is predicted, keep the results computed with that bit for the next round
                if bit == 1:
                    currentSet = set1[:]
                else:
                    currentSet = set0[:]

                if (error == 0) & (abs(diff1-diff2) > 50 ) :
                    stableSet = currentSet[:]

                # add the predicted bit to the key
                # print j+1, ": ",  bit, "------", abs(int(diff1 - diff2)), "--------", int(diff1), int(diff2), "error:", error
                d+= str(bit)
                (valid, bit ) = checkKey(mess_origin[0], results[0], d)
                if valid:
                    return d + str(bit)
        else :
            keySize +=16
            sampleSize +=1000
    return d




if ( __name__ == "__main__" ) :

    if len(sys.argv) < 3 :
      raise Exception("not enough argv")

    inputFile = open(sys.argv[2])
    n = int(inputFile.readline(), 16)
    e = int(inputFile.readline(), 16)
    inputFile.close()

  # Produce a sub-process representing the attack target.
    target = subprocess.Popen( args   = sys.argv[ 1 ],
                             stdout = subprocess.PIPE,
                             stdin  = subprocess.PIPE )

  # Construct handles to attack target standard input and output.
    target_out = target.stdout
    target_in  = target.stdin

    key = "1"

    global ro
    ro   = calc_ro(n)

    global ro_s
    ro_s = calc_ro_square(n)

    global omega
    omega = getOmega(n)
    # Execute a function representing the attacker.
    key = attack(n, e)
    print "FOUND KEY = ", "{0:X}".format(int(key, 2))
