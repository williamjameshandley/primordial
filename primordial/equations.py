class Equations(object):
    def sol(self, sol):
        return sol

class Event(object):
    def __init__(self, equations, direction=0, terminal=False):
        self.equations = equations
        self.direction = direction
        self.terminal = terminal
