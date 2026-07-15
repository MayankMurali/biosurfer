import sys
import traceback
from bisect import bisect
from collections.abc import Mapping
from contextlib import AbstractContextManager
from dataclasses import dataclass, fields
from operator import itemgetter
from typing import Generic, Iterable, Iterator, Tuple, TextIO, TypeVar

T = TypeVar('T')


class BisectDict(Mapping, Generic[T]):
    def __init__(self, items: Iterable[Tuple[int, T]]):
        self.breakpoints, self._values = zip(*sorted(items, key=itemgetter(0)))

    def __getitem__(self, key: int) -> T:
        if key < 0:
            raise KeyError('Key must be non-negative')
        i = bisect(self.breakpoints, key)
        try:
            return self._values[i]
        except IndexError as e:
            raise KeyError(key) from e

    def __iter__(self) -> Iterator[int]:
        yield from (0,) + self.breakpoints[:-1]

    def __len__(self):
        return len(self.breakpoints)


def frozendataclass(cls):
    frozencls = dataclass(cls, frozen=True)
    field_names = {field.name for field in fields(frozencls)}
    def replace(self, **kwargs):
        """Return new instance of frozendataclass with updated values."""
        new_field_values = {name: kwargs.get(name, getattr(self, name)) for name in field_names}
        return frozencls(**new_field_values)
    frozencls.replace = replace
    return frozencls


class ExceptionLogger(AbstractContextManager):
    def __init__(self, info=None, output: TextIO = None, callback=None):
        self.info = info
        self.callback = callback if callable(callback) else None
        self.output = output if output is not None else sys.stderr

    def __exit__(self, exc_type, exc_val, exc_tb):
        if exc_type is not None:
            self.output.write('---------\n')
            if self.info:
                self.output.write(str(self.info) + '\n')
            traceback.print_exc(file=self.output)
            self.output.write('---------\n')
            if self.callback:
                self.callback(exc_type, exc_val, exc_tb)
            return True
