from sage.crypto.boolean_function import BooleanFunction
from sage.crypto.sbox import SBox
import itertools
import random

os.environ['SBOX_ROOT'] = os.getcwd()
os.chdir(os.getcwd()+"/..")

load(os.environ['SBOX_ROOT'] + "/Cython/CFunc.spyx")
load(os.environ['SBOX_ROOT'] + "/Cython/CPPFunc.spyx")
load(os.environ['SBOX_ROOT'] + "/Sage/GSbox.sage")
load(os.environ['SBOX_ROOT'] + "/Sage/CSbox.sage")

n = 8
c = 65
gf = GF(2^(n-1))

binaryToInt = lambda binaryVector: int("".join(str(x) for x in binaryVector), 2)
getBin = lambda x, n: list(map(lambda x: int(x), list(format(x, 'b').zfill(n))))

def T(x1, x2, c, n):
	gf_x1 = gf.fetch_int(binaryToInt(x1))
	if gf_x1 == 0:
		return gf_x1
	reverse_x1 = gf_x1^(-1)
	return reverse_x1 * gf.fetch_int(c)^x2

def U(f, x2, x1):
	return (f(x1) + x2) % 2

def F(x):
	x1 = x[:-1]
	x2 = x[-1]
	T_res = getBin(T(x1, x2, c, n).integer_representation(), n-1)
	U_res = [U(randBoolFunc, x2, T_res)]
	#print("TRES: ", T_res)
	#print("URES: ", U_res)
	return T_res + U_res

randBoolFunc = BooleanFunction(list(map(lambda x: random.randint(0, 1), range(0,2^(n-1)))))

combinations = [list(i) for i in itertools.product([0, 1], repeat=n)]

def createTUSBoxes():
	cArray = []
	fixedPointsArray = []
	minDedreeArray = []
	sigmaArray = []
	trC_Array = []
	trCInv_Array = []
	degTArray = []
	degUArray = []
	AIArray = []

	for i in range(2,30):
		print("\nIteration: {}".format(i))
		c = i
		F_results = list(map(lambda x: F(x), combinations))
		T_results = list(map(lambda x: getBin(T(x[:-1], x[-1], c, n).integer_representation(), n-1), combinations))
		U_results = list(map(lambda x: U(randBoolFunc, x[-1], getBin(T(x[:-1], x[-1], c, n).integer_representation(), n-1)), combinations))

		s = list(map(lambda x: binaryToInt(x), F_results))
		T_s = map(lambda x: binaryToInt(x), T_results)
		U_func = BooleanFunction(U_results)

		sbox = Sbox(n=8,m=8,sbox=s.copy())
		sageBox = SBox(s.copy())
		T_box = SBox(T_s)

		cArray.append(c)
		print("c = ", c)

		fixedPoins = sbox.fixed_points()
		fixedPointsArray.append(fixedPoins)
		print("Fixed points: ", fixedPoins)

		#print("Sage Fixed points: ", sage_sbox.fixed_points())

		minDegree = sbox.minimum_degree()
		minDedreeArray.append(minDegree)
		print("Min degree: ", minDegree)

		#print("Polynimials lenght: ", len(sageBox.polynomials()))

		sigma = sageBox.differential_uniformity()
		sigmaArray.append(sigma)
		print("sigma: ", sigma)

		trC = gf.fetch_int(c).trace()
		trC_Array.append(trC)
		print("tr(c): ", trC)

		trCInc = (gf.fetch_int(c)^(-1)).trace()
		trCInv_Array.append(trCInc)
		print("tr(c^-1): ", trCInc)

		degT = T_box.min_degree()
		degTArray.append(degT)
		print("deg(T) = ", degT)

		degU = U_func.algebraic_degree()
		degUArray.append(degU)
		print("deg(U) = ", degU)

		ai = sbox.algebraic_immunity_sbox()
		AIArray.append(ai)
		print(("Algebraic immunity: degree={0} equations={1}".format(ai[0],ai[1])))

	print("\n\n ----------Results------------\n\n")
	print("rand f(x) degree: ", randBoolFunc.algebraic_degree(), "\n")
	print("cArray = ", cArray)
	print("fixedPointsArray = ", fixedPointsArray)
	print("minDedreeArray = ", minDedreeArray)
	print("sigmaArray = ", sigmaArray)
	print("trC_Array = ", trC_Array)
	print("trCInv_Array = ", trCInv_Array)
	print("degTArray = ", degTArray)
	print("degUArray = ", degUArray)
	print("AIArray = ", AIArray)
	print(randBoolFunc)
