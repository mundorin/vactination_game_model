class SIR:
    def __init__(self, s, i, r, b, y, n=10000, treshhold=0.0001):
        if s + i + r != 1:
            print("error s + i + r != 1")
        self.s = s
        self.i = i
        self.r = r
        self.b = b
        self.y = y
        self.time_step = 1
        self.timestamp = 0
        self.n = n
        self.treshhold = treshhold
        if self.i == 0:
            self.end = True
        else:
            self.end = False

    def iteretion(self):
        self.timestamp += self.time_step

        recovereds = self.y * self.i

        self.i -= recovereds
        self.r += recovereds

        sick = self.b * self.s * self.i

        self.s -= sick
        self.i += sick
        if not self.end and self.i <= self.treshhold:
            self.end = True

    def get_info(self):
        info = {
            "s": self.s,
            "i": self.i,
            "r": self.r,
            "b": self.b,
            "y": self.y,
            "n": self.n,
            "timestamp": self.timestamp,
            "end": self.end,
        }
        return info

class SIR_V:
    def __init__(self, s, i, r, v, b, y, n = 10000, treshhold = 0.0001):
        if s + i + r + v != 1:
            print("error s + i + r + v != 1")
        self.s = s
        self.i = i
        self.r = r
        self.v = v
        self.b = b
        self.y = y
        self.time_step = 1
        self.timestamp = 0
        self.n = n
        self.treshhold = treshhold
        if self.i == 0:
            self.end = True
        else:
            self.end = False

    def iteretion(self):
        self.timestamp += self.time_step

        recovereds = self.y * self.i

        self.i -= recovereds
        self.r += recovereds

        sikings = self.b * self.s * self.i

        self.s -= sikings
        self.i += sikings
        if not self.end and self.i <= self.treshhold:
            self.end = True


    def get_info(self):
        info = {
            "s": self.s,
            "i": self.i,
            "r": self.r,
            "b": self.b,
            "y": self.y,
            "v": self.v,
            "n": self.n,
            "timestamp": self.timestamp,
            "end": self.end,
        }
        return info

class SIRV:
    def __init__(self, s, i, r, v, b, y, fi, n = 10000, treshhold = 0.0001):
        if s + i + r + v != 1:
            print("error s + i + r + v != 1")
        self.s = s
        self.i = i
        self.r = r
        self.v = v
        self.b = b
        self.y = y
        self.fi = fi
        self.time_step = 1
        self.timestamp = 0
        self.n = n
        self.treshhold = treshhold
        if self.i == 0:
            self.end = True
        else:
            self.end = False

    def iteretion(self):
        self.timestamp += self.time_step

        vactinatons = self.s * self.fi

        self.s -= vactinatons
        self.v += vactinatons

        recovereds = self.y * self.i

        self.i -= recovereds
        self.r += recovereds

        sikings = self.b * self.s * self.i

        self.s -= sikings
        self.i += sikings

        if not self.end and self.i <= self.treshhold:
            self.end = True


    def get_info(self):
        info = {
            "s": self.s,
            "i": self.i,
            "r": self.r,
            "b": self.b,
            "y": self.y,
            "v": self.v,
            "fi": self.fi,
            "n": self.n,
            "timestamp": self.timestamp,
            "end": self.end,
        }
        return info







