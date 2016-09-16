import sys, subprocess, hashlib, shlex, math, binascii, random, copy
import threading
from multiprocessing import Process


global ro
global ro_s
global omega

global sampleSize
sampleSize = 1500

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

def calc_ro(nrOfBits):
    ro = 1
    while ro < nrOfBits:
        ro = ro << 64
    return  ro

def interact(  c ) :
  global interactions
  interactions += 1
  target_in.write("{0:X}".format(c) + "\n")

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

    mont_messages = []
    messages = []
    timings = []
    results = []
    for j in range(0, sampleSize):
        m = random.getrandbits(n.bit_length())
        while m > n :
            m = random.getrandbits(n.bit_length())
        messages.append(m)
        mont_messages.append(MontMul(messages[j], ro_s, n, omega)[0])

        res = interact(m)
        timings.append(res[0])
        results.append(res[1])
    return (messages, mont_messages, timings, results)

def updateMessages(messages, currentKey, n):
    currentSet = messages[:]
    for i in range(0, len(messages)):
        currentSet[i] = MontMul(currentSet[i], currentSet[i], n, omega)[0]

    for i in range(0, len(messages)):
        for j in range(1, len(currentKey)):
            if(currentKey[j] == "1"):
                currentSet[i] = MontMul(currentSet[i], messages[i], n, omega)[0]
            currentSet[i] = MontMul(currentSet[i], currentSet[i], n , omega)[0]
    return currentSet

def checkKey( message, result, d):
    key_1 = int(d + "1", 2)
    key_0 = int(d + "0", 2)

    # choose which prediction was right and output key
    if pow(message, key_1, nrOfBits) == result:
        return (1, "1")

    if pow(message, key_0, nrOfBits)== result:
        return (1, "0")

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
    def __init__(self, nrOfBits, sampleSize):
        threading.Thread.__init__(self)
        self.n = nrOfBits
        self.sampleSize = sampleSize
    def run(self):
        ro_s = calc_ro_square(self.n)
        omega = getOmega(self.n)
        messages = []
        while len(messages) < self.sampleSize:
            m = random.getrandbits(nrOfBits.bit_length())
            while m > nrOfBits :
                m = random.getrandbits(nrOfBits.bit_length())
            messages.append(m)
        self.messages = messages
    def join(self):
        return self.messages

def checkAnormal(diff1, diff2):
    if (abs(diff1 - diff2) < treshold) | ((diff1 < 0) & (diff2 < 0)):
        return 1
    return 0

class Group(object):
    time = 0
    size = 0

def simulateStep( mont_messages, currentSet, nrOfBits, groups, time  ):
    encoded_1 = currentSet[:]
    encoded_0 = currentSet[:]
    for i in range(0, len(mont_messages)):
        # assume bit j is 0
        (encoded_0[i], extra) = MontMul(currentSet[i], currentSet[i], nrOfBits, omega)
        # if extra reduction
        if extra :
            groups[3].time += time[i]
            groups[3].size += 1
        # if no extra reduction
        else :
            groups[4].time += time[i]
            groups[4].size += 1

        # assume bit j is 1
        (temp, extra) = MontMul(currentSet[i], mont_messages[i], nrOfBits , omega)
        (encoded_1[i], extra) = MontMul(temp, temp, nrOfBits ,omega)

        # if extra reduction
        if extra :
            groups[1].time += time[i]
            groups[1].size +=1
        # if no extra reduction
        else :
            groups[2].time += time[i]
            groups[2].size +=1
    return (groups, encoded_0, encoded_1)

def compute_differences( groups ):

    # compute averaget time for all 4 groups
    uF1 = float(groups[1].time)/groups[1].size
    uF2 = float(groups[2].time)/groups[2].size
    uF3 = float(groups[3].time)/groups[3].size
    uF4 = float(groups[4].time)/groups[4].size

    # compute differences between pari groups
    diff1 = uF1 - uF2
    diff2 = uF3 - uF4

    return (diff1, diff2)

def look_Ahead(bitCheck, encoded_0, encoded_1, messages, results, mont_messages, nrOfBits, time, key):
    global lookAhead

    key += str(bitCheck)
    (validKey, bit) = checkKey(messages[0], results[0], key)
    if validKey :
        key += str(bit)
        return (validKey, key, delta)

    if bitCheck :
        temp_encoded  = encoded_1[:]
    else:
        temp_encoded  = encoded_0[:]
    delta = 0
    for y in range(0, lookAhead):
        groups  = [ Group() for i in range(5)]

        (groups, encoded_0, encoded_1 ) = simulateStep( mont_messages, temp_encoded, nrOfBits, groups, time );

        (diff1, diff2) = compute_differences( groups )

        bit = 0
        if diff1 > diff2:
            bit = 1

        if bit :
            temp_encoded = encoded_1[:]
        else:
            temp_encoded = encoded_0[:]
        key += str(bit)
        # print "ver1 : ", len(key_0), ": ",  bit, "------", abs(int(diff1 - diff2)), "--------", int(diff1), int(diff2), "error:", error
        (validKey, bit) = checkKey(messages[0], results[0], key)
        if validKey :
            key += str(bit)
            return (validKey, key, delta)

        if not checkAnormal(diff1, diff2) :
            delta += abs(int(diff1 - diff2))

    return (validKey, key, delta)

def attack(nrOfBits, e):
    global sampleSize

    # set start conditions
    validKey = 0
    keySize  = 64
    foundKey = "1"

    # while the key is within the limits of the maximum key and it is not valid
    while (keySize < maxKeySize) &  (not validKey):

        # re-initialise everything
        messages        = []
        mont_messages   = []
        currentSet      = []
        results         = []
        time            = []
        stablekeySet    = 0
        error           = 0

        # generate new sample messages of nrOfBits bits
        generateMessages = sampleMessages(nrOfBits, sampleSize)
        generateMessages.run()

        # add new messages to the current set of messages
        messages += generateMessages.join()

        # reset key to stable key
        currentKey = foundKey

        # compute times taken to decrypt each message
        for i in range(0, len(messages)):

            # get times for each message ( res[0] is the time taken to decrypt, and res[1] is decrypted message)
            res = interact(messages[i])
            results.append(res[1])
            time.append( res[0] )

            # transform each message into montgomery form
            mont_messages.append(MontMul(messages[i], ro_s, nrOfBits, omega)[0])

        print "Sample Size: ", sampleSize
        print "Interactions: ", interactions
        print "Start from key bit ", len(currentKey)

        currentSet = updateMessages(mont_messages, currentKey, nrOfBits)
        encoded_1 = currentSet[:]
        encoded_0 = currentSet[:]

        # discover the next bits untill the last bit which needs to be guessed
        while ((len( currentKey ) <= keySize) & (error < 7)):
            warning = 0
            groups  = [ Group() for i in range(5)]

            (groups, encoded_0, encoded_1 ) = simulateStep( mont_messages, currentSet, nrOfBits, groups, time );

            (diff1, diff2) = compute_differences( groups )

            # condition for anormal behaviour
            if checkAnormal(diff1, diff2) :
                warning = 1
                error += 2
            else :
                if error > 0:
                    error -= 1
                warning = 0

            if warning > 0 :
                print "warning at bit", len( currentKey )
                stablekeySet = 1

                # check following rounds for bit 0
                (validKey, possibleKey, delta0) = look_Ahead(0, encoded_0, encoded_1, messages, results, mont_messages, nrOfBits, time, currentKey)
                if(validKey):
                    return possibleKey

                # check following rounds for bit 1
                (validKey, possibleKey, delta1) = look_Ahead(1, encoded_0, encoded_1, messages, results, mont_messages, nrOfBits, time, currentKey)
                if(validKey):
                    return possibleKey

                # decide which set was better
                if(delta0 < delta1) :
                    bit = 1
                else :
                    bit = 0

            else :
                #if diff for the groups with assumed bit = 1 is bigger than the diff of the groups with assumed bit = 0, then predict 1, else predict 0
                if diff1 > diff2:
                    bit = 1
                else :
                    bit = 0

                if(not stablekeySet):
                    foundKey = currentKey

            # depending on which bit is predicted, keep the results computed with that bit for the next round
            if bit == 1:
                currentSet = encoded_1[:]
            else:
                currentSet = encoded_0[:]

            # add the predicted bit to the key
            currentKey += str(bit)
            (validKey, bit ) = checkKey(messages[0], results[0], currentKey)
            if validKey:
                return currentKey + str(bit)

        else :
            keySize += 16
            sampleSize += 1000
    return currentKey

if ( __name__ == "__main__" ) :

    if len(sys.argv) < 3 :
      raise Exception("not enough argv")

    inputFile = open(sys.argv[2])
    nrOfBits = int(inputFile.readline(), 16)
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
    ro   = calc_ro(nrOfBits)

    global ro_s
    ro_s = calc_ro_square(nrOfBits)

    global omega
    omega = getOmega(nrOfBits)
    # Execute a function representing the attacker.
    key = attack(nrOfBits, e)
    print "FOUND KEY = ", "{0:X}".format(int(key, 2))
