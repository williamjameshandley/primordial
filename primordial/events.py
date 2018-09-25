class Event(object):
    def __init__(self, equations, direction=0, terminal=False):
        self.equations = equations
        self.direction = direction
        self.terminal = terminal


class Inflation(Event):
    def __call__(self, t, y):
        return self.equations.inflating(y)


class Stationary(Event):
    def __call__(self, t, y):
        return self.equations.H2(y)

