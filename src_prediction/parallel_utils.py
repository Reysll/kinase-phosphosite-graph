from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from typing import Callable, Iterable, List, TypeVar


T = TypeVar("T")
R = TypeVar("R")


def run_in_parallel(
    items: Iterable[T],
    worker_fn: Callable[[T], R],
    max_workers: int,
    use_threads: bool = False,
) -> List[R]:
    """
    Run worker_fn over items in parallel and return results in completion order.

    use_threads=True  — ThreadPoolExecutor: threads share memory (no copying of
                        large numpy arrays). sklearn/liblinear releases the GIL
                        during fitting, so real speedup is achievable. Best for
                        single-machine runs where large arrays are pre-built.

    use_threads=False — ProcessPoolExecutor: separate processes. Required when
                        worker_fn uses code that does not release the GIL.
                        On Linux (fork), pre-built arrays are copy-on-write shared
                        so memory duplication is minimal. On Windows (spawn), large
                        arrays in tasks are pickled per task — avoid for big data.
    """
    items = list(items)
    if max_workers <= 1:
        return [worker_fn(item) for item in items]

    executor_cls = ThreadPoolExecutor if use_threads else ProcessPoolExecutor
    results: List[R] = []
    with executor_cls(max_workers=max_workers) as ex:
        futures = {ex.submit(worker_fn, item): item for item in items}
        for fut in as_completed(futures):
            results.append(fut.result())
    return results
