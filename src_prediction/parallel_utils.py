from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from typing import Callable, Iterable, List, Optional, TypeVar


T = TypeVar("T")
R = TypeVar("R")


def run_in_parallel(
    items: Iterable[T],
    worker_fn: Callable[[T], R],
    max_workers: int,
    use_threads: bool = False,
    progress_fn: Optional[Callable[[int, int], None]] = None,
) -> List[R]:
    """
    Run worker_fn over items in parallel and return results in completion order.

    progress_fn(done, total) is called after each task completes, from inside
    the as_completed loop, so timing/ETA reflect actual wall-clock progress.

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
    n = len(items)

    if max_workers <= 1:
        results: List[R] = []
        for i, item in enumerate(items, 1):
            results.append(worker_fn(item))
            if progress_fn:
                progress_fn(i, n)
        return results

    executor_cls = ThreadPoolExecutor if use_threads else ProcessPoolExecutor
    results = []
    with executor_cls(max_workers=max_workers) as ex:
        futures = {ex.submit(worker_fn, item): item for item in items}
        for i, fut in enumerate(as_completed(futures), 1):
            results.append(fut.result())
            if progress_fn:
                progress_fn(i, n)
    return results
