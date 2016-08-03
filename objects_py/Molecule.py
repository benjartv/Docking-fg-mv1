__author__ = 'Felipe_Gonzalez_Foncea'

from utils_py.BiopyVector import *
from utils_py.Methods import *


#------------------------------ Molecule ----------------------------
class Molecule(object):


    #-------------------------------- Init ------------------------------
    def __init__(self, recordName = ""):
        self.recordName = recordName
        self.atom = []
        self.alterLocation = []
        self.residueName = []
        self.chainId = []
        self.residueSeq = []
        self.insertCode = []
        self.x = []
        self.y = []
        self.z = []
        self.occupancy = []
        self.tempFactor = []
        self.segmentId = []
        self.symbol = []
        self.charge = []
        self.connectMatrix = []
        self.rotationsPoints = []
        self.rotationsSegment = []


    #------------------------------ Traslate ----------------------------
    def translate(self, delta):
        for i in xrange(len(self.x)):
            self.x[i] += delta[0]
            self.y[i] += delta[1]
            self.z[i] += delta[2]


    #-------------------------- Traslate to Point -----------------------
    def traslateToPoint(self, point):
        originPoint = self.findCenter()
        delta = [point[0] - originPoint[0],
                 point[1] - originPoint[1],
                 point[2] - originPoint[2]]

        for i in xrange(len(self.x)):
            self.x[i] += delta[0]
            self.y[i] += delta[1]
            self.z[i] += delta[2]


    #------------------------- Add Reference Point ----------------------
    def referencePoint(self, x, y, z, namePoint):
        self.atom.append("O")
        self.alterLocation.append("")
        self.residueName.append(namePoint)
        self.chainId.append("A")
        self.residueSeq.append("999")
        self.insertCode.append("")
        self.x.append(x)
        self.y.append(y)
        self.z.append(z)
        self.occupancy.append("1.00")
        self.tempFactor.append("0.0")
        self.segmentId.append("")
        self.symbol.append("O")
        self.charge.append("")


    #-------------------------- Rotate by Angle -------------------------
    def rotateByAngle(self, pointA, pointB, theta):
        for j in xrange(len(self.x)):
            self.x[j] -= pointA[0]
            self.y[j] -= pointA[1]
            self.z[j] -= pointA[2]

        p1 = [pointA[0], pointA[1], pointA[2]]
        p2 = [pointB[0], pointB[1], pointB[2]]

        vetorRef = Vector(p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2])
        matrixRot = rotaxis2m(theta, vetorRef)

        for i in xrange(len(self.x)):
            pointRef = [[self.x[i]],[self.y[i]],[self.z[i]]]
            mRotacionada = crossProduct(matrixRot, pointRef)
            self.x[i] = mRotacionada[0][0] + pointA[0]
            self.y[i] = mRotacionada[1][0] + pointA[1]
            self.z[i] = mRotacionada[2][0] + pointA[2]


     #--------------------------- Rotate by Vector -----------------------
    def rotateByVector(self, vector, theta):
        vector = self.validateNormCero(vector)

        originPoint = self.findCenter()

        unitX = vector[0] + originPoint[0]
        unitY = vector[1] + originPoint[1]
        unitZ = vector[2] + originPoint[2]

        for j in xrange(len(self.x)):
            self.x[j] -= originPoint[0]
            self.y[j] -= originPoint[1]
            self.z[j] -= originPoint[2]

        p1 = [originPoint[0], originPoint[1], originPoint[2]]
        p2 = [unitX, unitY, unitZ]

        vetorRef = Vector(p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2])
        matrixRot = rotaxis2m(theta, vetorRef)

        for i in xrange(len(self.x)):
            pointRef = [[self.x[i]],[self.y[i]],[self.z[i]]]
            mRotacionada = crossProduct(matrixRot, pointRef)
            self.x[i] = mRotacionada[0][0] + originPoint[0]
            self.y[i] = mRotacionada[1][0] + originPoint[1]
            self.z[i] = mRotacionada[2][0] + originPoint[2]


    #--------------------------- Rotate Segment -------------------------
    def rotateSegment(self, idRotation, theta):
        serialA = self.rotationsPoints[idRotation][0]
        serialB = self.rotationsPoints[idRotation][1]
        delthaX = self.x[serialA]
        delthaY = self.y[serialA]
        delthaZ = self.z[serialA]

        for j in xrange(len(self.x)):
            self.x[j] -= delthaX
            self.y[j] -= delthaY
            self.z[j] -= delthaZ

        p1 = [self.x[serialA], self.y[serialA], self.z[serialA]]
        p2 = [self.x[serialB], self.y[serialB], self.z[serialB]]

        vetorRef = Vector(p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2])
        matrixRot = rotaxis2m(theta, vetorRef)

        segment = []
        segment.append(serialA)
        segment = self.findRotationsSegment(segment, self.connectMatrix,serialB)
        self.rotationsSegment = segment[:]

        for i in self.rotationsSegment:
            pointRef = [[self.x[i]], [self.y[i]], [self.z[i]]]
            mRotacionada = crossProduct(matrixRot, pointRef)
            self.x[i] = mRotacionada[0][0]
            self.y[i] = mRotacionada[1][0]
            self.z[i] = mRotacionada[2][0]

        for j in xrange(len(self.x)):
            self.x[j] += delthaX
            self.y[j] += delthaY
            self.z[j] += delthaZ


    #------------------------ Find Rotation Segment ---------------------
    def findRotationsSegment(self, result, matrix, id):
        for i in xrange(len(matrix[id])):
            if(matrix[id][i]==1 and i not in result ):
                result.append(i)
                aux = self.findRotationsSegment(result, matrix, i)
                if(len(aux) <> 0):
                    result = aux[:]
        return result[:]


    #----------------------- Find center of Molecule --------------------
    def findCenter(self):
        sumX, sumY, sumZ = 0, 0, 0
        center = [0,0,0]
        sizeM = len(self.x)

        if (sizeM > 0):
            for i in xrange(sizeM):
                sumX += self.x[i]
                sumY += self.y[i]
                sumZ += self.z[i]
            center[0] = sumX / sizeM
            center[1] = sumY / sizeM
            center[2] = sumZ / sizeM

        return center


    #-------------------- Find center atom of Molecule -------------------
    def findCenterAtom(self):
        sumX, sumY, sumZ = 0, 0, 0
        center = [0,0,0]
        sizeM = len(self.x)
        minDistance = 100
        lastPoint = 0

        if (sizeM > 0):
            for i in xrange(sizeM):
                sumX += self.x[i]
                sumY += self.y[i]
                sumZ += self.z[i]

            center[0] = sumX / sizeM
            center[1] = sumY / sizeM
            center[2] = sumZ / sizeM

            for i in xrange(sizeM):
                point = [self.x[i], self.y[i], self.z[i]]
                distance = euclideanDistance(center,point)

                if (distance < minDistance):
                    minDistance = distance
                    lastPoint = point
        return lastPoint


    #------------------------- Find Serial Number -----------------------
    def findSerialNumber(self, atom):
        if (type(atom) is str):
            result = -1
            for i in xrange(len(self.x)):
                if(atom.strip() == self.atom[i].strip()):
                    result = i
                    break
        else:
            result = atom - 1
        return result


    #------------------------- Add Rotation Atoms -----------------------
    def addRotationsAtoms(self, atomsReferences):
        for ref in atomsReferences:
            self.rotationsPoints.append([self.findSerialNumber(ref[0]), self.findSerialNumber(ref[1])])


    #------------------------ Validate Norm Cero ------------------------
    def validateNormCero(self,vector):
        res = 0.0
        res += vector[0] * vector[0]
        res += vector[1] * vector[1]
        res += vector[2] * vector[2]
        res = math.sqrt(res)

        if(res<0.2):
            vector = randomSpherePoint(1)
        return vector
