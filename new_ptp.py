import numpy as np
from numpy import linalg as LA
from choldate import cholupdate
import time
from scipy import linalg
from decimal import Decimal

times = []

def Timer(f): # Функция таймер, в которую будет "обернута" будущая функция
    def func(*args, **kwargs): # Функция которая будет выполнять переданную таймеру функцию и засекать время, требуется для того
        # чтоб выполнялась функция с произвольными начальными данными
        currtime = time.time() # Запомнить текущее время
        result = f(*args, **kwargs) #Выполнить введенную в таймер функцию
        times.append(time.time()-currtime)
        print ("Time:", (time.time()-currtime)) # Вывести результат
        return result

    return func # Вернуть функцию

#@profile 
def cake (dim, nvec, scale, shift) :
    np.random.seed(dim+nvec)
    X = np.random.uniform(0, 1, (dim - 1, nvec))
    lower_row = scale[1]*np.random.uniform(0,1, (1, nvec)) + shift
    vector = np.sum(X, axis=1)/nvec
    upper = np.dot( np.reshape(vector, (vector.shape[0], 1)) , np.ones((1, nvec)) )
    X = np.r_[ scale[0]*(X - upper), lower_row ]
    return X

#@profile 
def  cIn(r, k, u) :

    c = linalg.solve_triangular(r.T, u[:get_row(u)-1,:] , lower=True)

    if ( u[get_row(u)-1,0] - np.dot(c.T,c) )<0 :
        print("PANIC:", ( u[get_row(u)-1,0] - np.dot(c.T,c) ) )
        raise Exception('Incorrect data')

    d = np.sqrt( u[get_row(u)-1,0] - np.dot(c.T,c) )

    newR = np.c_[ r, c ]
    down = np.zeros((1, newR.shape[1]))
    down[0, get_column(down)-1] = d
    newR = np.r_[ newR, down ]

    return newR

#@profile
def  newDel(r, iout) :

    lowMinor = np.copy(r)[iout+1:, iout+1:]
    row = np.copy(r)[iout,:]
    row = row[iout+1:]

    cholupdate(lowMinor, row)

    newR = np.copy(r)
    newR[iout+1:, iout+1:] = lowMinor
    newR = np.delete(newR, (iout), axis=0)
    newR = np.delete(newR, (iout), axis=1)

    return newR

#@profile 
def sumsq(x):
    sq = np.zeros((x.shape[1]))
    for index, value in enumerate(x.T):
        sq[index] = np.inner(value, value)
    return sq
#@profile 
def get_min(a) :
    if a.ndim == 2 and a.shape[0]==1 :
        a = np.reshape(a, a.shape[1],0)
    minimum = np.amin(a, axis=0)
    index = np.argmin(a, axis=0)
    return minimum, index
#@profile 
def first_iteration(vectors):

    vmin, ifac = get_min(sumsq(vectors))
    kvec = np.array([ifac])
    lmb = np.array([ 1 ])
    r = LA.norm(vectors[:, ifac])
    z = vectors[:, kvec]*lmb

    return vmin, ifac, kvec, lmb, r, z

def kvec_iteration(vectors, kvec):

    lmb = np.ones((len(kvec), 1))*(1/len(kvec))
    r = np.linalg.cholesky(np.dot(vectors[:, kvec].T, vectors[:, kvec]))
    r = r.T
    z = np.dot(vectors[:, kvec], lmb)
    return lmb, r, z


#@profile
def getNewVector(vectors, z, accuracy) :

    stop = False
    get_ifac_var = np.dot(z.T, vectors) - sumsq(z)
    vmin, ifac = get_min(get_ifac_var)
    reps = accuracy*LA.norm(vectors[:, ifac])

    if ( vmin > -reps ).all() :
        stop = True

    return vmin, ifac, stop

#@profile
def report(addIter, remIter, add, remove, z, vmin, kvecLen, iteration_tick) :

    if len(remove) == 0:
        remove = ["-"]
    '''
    template = "{0:>7}{1:6}{2:>3}{3:6}{4:>3}{5:>4}{6:6}{7:>4}{8:>15}{9:>4}{10:>6}{11:>4}{12:>6}{13:>4}{14:>15}{15:>5}{16:>10}"
    print (template.format("++ iter", addIter,
                           "(+)", remIter,
                           "(-)", "len",
                           kvecLen, "zx",
                           format(vmin, '.4f'), "in",
                           " ".join(str(i) for i in add), "out",
                           " ".join(str(i) for i in remove),  "zz",
                           format(sumsq(z)[0], '.4f'), "tic",
                           format(time.time()-iteration_tick, '.4f')             ) )
    '''

    template = "{0:>7}{1:6}{2:>3}{3:6}{4:>3}{5:>4}{6:6}{7:>4}{8:>12.4E}{9:>4}{10:>6}{11:>4}{12:>6}{13:>4}{14:>12.4E}{15:>5}{16:>8}"
    print (template.format("++ iter", addIter,
                       "(+)", remIter,
                       "(-)", "len",
                       kvecLen, "zx",
                       Decimal(vmin), "in",
                       " ".join(str(i) for i in add), "out",
                       " ".join(str(i) for i in remove),  "zz",
                       Decimal(sumsq(z)[0]), "tic",
                       format(time.time()-iteration_tick, '.4f')             ) )


    iteration_tick = time.time()

    #print(addIter+remIter, "+", addIter, "-", remIter, "Add", add, "Rem", remove, "vmin", format(vmin, '.4f'), "zz", format(sumsq(z)[0], '.4f'))
    return iteration_tick

#@profile 
def get_column(a):
    if a.ndim == 1:
        return a.shape[0]
    else:
        return a.shape[1]
#@profile 
def get_row(a):
    if a.ndim == 1:
        return 1
    else:
        return a.shape[0]

#@profile 
def addVector (kvec, lmb, ifac) :

    kvec = np.append(kvec, ifac)
    lmb = np.append(lmb, 0)
    return kvec, lmb

#@profile 
def calculateChol (X, kvec, r, ifac) :

    if get_column(kvec) <= get_row(X) :

        u = np.dot(X[:, kvec].T, X[:, ifac])
        u = np.reshape(u, (-1, 1))

        newR = cIn(r, get_column(kvec), u)
        return newR

    return r

#@profile 
def project(X, kvec, r, z) :

    if get_column(kvec) <= get_row(X) :
        lmb = linalg.solve_triangular( r ,   linalg.solve_triangular ( (r).T , np.ones((get_column(kvec), 1)), lower=True ) )
        lmb /= sum(lmb)
        z = np.dot(X[:, kvec], lmb)
    else:
        xi = X[:, kvec[len(kvec)-1]]
        xi = np.reshape(xi, (-1, 1))
        mu_last = sumsq(z)/(sumsq(z) - np.dot(xi.T, z))
        xp = (mu_last/(1 - mu_last))*np.dot(X[:,kvec[0:len(kvec)-1]].T, xi)
        q = linalg.solve_triangular(r, linalg.solve_triangular(r.T, np.ones(z.shape), lower=True) )
        tt = -(1 + np.dot(q.T, xp))/np.dot(q.T, np.ones(z.shape))
        mu = tt + xp
        mu = linalg.solve_triangular(-r, linalg.solve_triangular(r.T, mu, lower=True))
        mu *= (1 - mu_last)
        mu = np.append(mu, mu_last)
        lmb = np.reshape(mu, (-1, 1))
        z = np.zeros(z.shape)

    return lmb, z

#@profile 
def bridge(vectors, r, kvec, lmb, z, mayBeLmb):

    indx = np.array(mayBeLmb < 0)
    indx = indx.astype(bool)
    if all(v == 0 for v in indx):
        raise ValueError("99")

    art = np.reshape(lmb, indx.shape)[indx] / (np.reshape(lmb, indx.shape) - mayBeLmb)[indx]

    theta, imin = get_min(art)
    iout = np.arange(get_row(mayBeLmb))[indx.ravel()]
    iout = iout[imin]

    auxlmb = theta*mayBeLmb + (1 - theta)*np.reshape(lmb, mayBeLmb.shape)

    new_lmb = np.append( auxlmb[0:iout], auxlmb[iout+1:len(auxlmb)] )
    new_lmb = np.reshape(new_lmb, (-1, 1))
    new_kvec = np.append(kvec[0:iout], kvec[iout+1:len(kvec)] )
    new_R = newDel(r, iout)

    if get_column(kvec) > get_row(X):
        ilast = kvec[len(kvec)-1]
        u = np.dot( X[:, new_kvec[0:len(new_kvec)-1]].T , X[:, ilast])
        u = np.append(u, sumsq(X[:, [ilast]]) )
        u = np.reshape(u, (-1, 1))
        attR = cIn(new_R, get_column(new_kvec), u)
        new_R = attR

    new_z = np.dot ( X[:, new_kvec], np.reshape(new_lmb,  (-1, 1)) )

    return new_kvec, new_lmb, new_z, new_R, iout

#@profile 
@Timer
def ptp (vectors, kvec, maximum_iteration_number, accuracy, pub) :

    iteration_tick = time.time()

    vectors = np.array(vectors)

#-----------first iteration
    addIter = 1
    remIter = 0
    stop = False


    if len(kvec)==0:
        vmin, ifac, kvec, lmb, r, z = first_iteration(vectors)

        if r == 0:
            raise  Exception('0 column in X')

        r = [[r]]
        r = np.array(r)

        if ((addIter+remIter)%pub == 0) :
            iteration_tick = report(addIter, remIter, kvec, [], z, vmin, len(kvec), iteration_tick)

    else:
        lmb, r, z = kvec_iteration(vectors, kvec)

    #-----------main body

    old_kvec = kvec

    while ( (not stop) and (addIter+remIter <= maximum_iteration_number) ) :

#add vector--------------------------------------

        vmin, ifac, stop = getNewVector(vectors, z, accuracy)

        if stop :
            print("__________________________________________________________________________________________________")
            iteration_tick = report(addIter, remIter, add, rem, z, vmin, len(kvec), time.time())
            return kvec, z, lmb
        else:
            kvec, lmb = addVector(kvec, lmb, ifac)
            r = calculateChol (vectors, kvec, r, ifac)

#project--------------------------------------

        mayBeLmb, mayBeZ = project(vectors, kvec, r, z)
        flag = False

#check project---------------------------------

        while any(mayBeLmb < 0) :
            kvec, lmb, z, r, iout = bridge(vectors, r, kvec, lmb, z, mayBeLmb)
            flag = True
            mayBeLmb, mayBeZ = project(vectors, kvec, r, z)

#true projection--------------------------------

        lmb, z = mayBeLmb, mayBeZ

        if flag:
            remIter+=1
        else:
            addIter+=1

        add = list(set(kvec) - set(old_kvec))
        rem = list(set(old_kvec) - set(kvec))
        old_kvec = kvec

        if ((addIter+remIter)%pub == 0) :
            iteration_tick = report(addIter, remIter, add, rem, z, vmin, len(kvec), iteration_tick)

    return kvec, z, lmb

def checkOptimal (vectors, z, accuracy):

    print ("Min check val:", get_min(np.dot(z.T,vectors)-sumsq(z)))

    if (np.dot(z.T,vectors)-sumsq(z) > -accuracy).all():
        return True
    else:
        return False


def getKvec (vectors, dim, vecNum):

    if dim>vecNum:
        num = vecNum-5
    else:
        num = dim -10

    if num <= 0 or num >= dim:
        return []
    kvec = sumsq(vectors).argsort()[:num]
    print("Initial Len:", len(kvec))
    print("Initial kvec:", kvec)

    return kvec

def secondGetKvec (vectors, dim, vecNum):

    if dim>vecNum:
        num = vecNum - 2
    else:
        num = dim - 2

    if num <= 0 or num >= dim:
        return []

    zf = np.sum(vectors, axis=1)
    zf = np.reshape(zf, (1, zf.shape[0]))

    kvec = np.array([])

    for i in range(0, vecNum):
        if np.dot(zf, vectors[:,i])<=LA.norm(zf):
            kvec = np.append(kvec, i)

    kvec = kvec.astype(int)
    return kvec [: num]

dim = 5000

for vecNum in range(5100, 5101, 50):
    print(dim,vecNum)
    X = 10 * cake(dim, vecNum, [50.0, 0.01], 0.0001)

    #print(X.shape)
    kvec = []
    #kvec = getKvec(X, dim, vecNum)
    #kvec = secondGetKvec(X, dim, vecNum)
    print("Kvec:", kvec)

    kvec, z, lmb = ptp(X, kvec, 10000, 10**(-10), 1)

    #print("###########################")
    #print("Kvec:", kvec)
    #print("Z:", z)
    #print("Lmb:", lmb)

    print (checkOptimal(X, z, 0.00001))
    print("###########################")

print (times)