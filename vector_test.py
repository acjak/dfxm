#!/bin/python
# -*- coding: utf-8 -*-
"""blah."""

from math import sin, cos, tan, sqrt

om = 80.
th = 10.

ax = sin(om) * (sin(th) / tan(th))
ay = -cos(om) * (sin(th) / tan(th))
az = -sin(th)

print ax, ay, az

print sqrt(ax**2 + ay**2 + az**2)

print sqrt(sin(om)**2 + cos(om)**2)
