__author__ = 'Felipe_Gonzalez_Foncea'

import datetime
import time

class TimeManager(object):

    startTime = time.time()
    stopTime = time.time()

    startDt = datetime.datetime.now()
    stopDt = datetime.datetime.now()

    totalTime = 0.

    def __init__(self):
        self.totalTime = 0.


    def start(self):
        self.startTime = time.time()
        self.startDt = datetime.datetime.now()


    def stop(self):
        self.stopTime = time.time()
        self.stopDt = datetime.datetime.now()
        self.totalTime = (self.stopTime - self.startTime)


    def clear(self):
        self.totalTime = 0.


    def getTimeSeconds(self):
        return self.totalTime


    def getTimeMinutes(self):
        return self.totalTime / 60


    def getTimeHours(self):
        return self.totalTime / 3600


    def getTimeChronometer(self):
        return (self.stopDt - self.startDt)


    def printTimeSeconds(self, text=''):
        print text + ("{0:.4f}".format(self.totalTime)) + " s"


    def printTimeMinuts(self, text=''):
        if (self.totalTime == 0):
            print text + ("{0:.4f}".format(0))
        else:
            print text + ("{0:.4f}".format(self.totalTime/60)) + " m"


    def printTimeHours(self, text=''):
        if (self.totalTime == 0):
            print text + ("{0:.4f}".format(0))
        else:
            print text + ("{0:.4f}".format(self.totalTime/3600)) + " h"


    def printTimeChronometer(self, text=''):
        if (self.totalTime == 0):
            print text + "00:00:00"
        else:
            print text + (self.stopDt - self.startDt)