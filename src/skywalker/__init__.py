"""Skywalker: plan and visualize altitude/azimuth tracks of astronomical
objects for observers at any observatory on Earth."""

from .cli import Skywalker, parse_args, main

__all__ = ["Skywalker", "parse_args", "main"]

__version__ = "0.1.0"
