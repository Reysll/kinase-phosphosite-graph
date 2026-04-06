from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Callable, Iterable, List, TypeVar


T = TypeVar("T")
R = TypeVar("R")


def run_in_parallel(
    items: Iterable[T],
    worker_fn: Callable[[T], R],
    max_workers: int,
) -> List[R]:
    """
    Run worker_fn over items in parallel processes and return results in completion order.
    """
    items = list(items)
    if max_workers <= 1:
        return [worker_fn(item) for item in items]

    results: List[R] = []
    with ProcessPoolExecutor(max_workers=max_workers) as ex:
        futures = {ex.submit(worker_fn, item): item for item in items}
        for fut in as_completed(futures):
            results.append(fut.result())
    return results