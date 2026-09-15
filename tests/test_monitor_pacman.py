"""Behavioral regressions for the decorative SIMU board."""

from jarvishep2.monitor.bars import visible_width
from jarvishep2.monitor.overview import PacmanGame


def board(width=12, height=2):
    game = PacmanGame()
    game._ensure_board(width, height)
    game._blocks = [set() for _ in range(height)]
    game._giant_spawn_ticks = 10000
    return game


def test_collision_uses_two_occupied_cells_and_checks_arrival():
    game = board()
    ghost = game._ghosts[0]
    game._ghosts = [ghost]
    ghost.row, ghost.column, ghost.speed = 0, 4, 2
    assert not game._handle_collisions(0, 3)
    game._path_index = 3
    game._step = 1  # This ghost stays still; Pac-Man moves into it.
    game._advance_game()
    assert game._path[game._path_index] == (0, 4)
    assert game._caught_ticks == game._CAUGHT_TICKS
    assert not game._pellets[0][4]


def test_row_end_catch_does_not_start_recovery_until_departure():
    game = board()
    ghost = game._ghosts[0]
    game._ghosts = [ghost]
    ghost.row, ghost.column = 0, 10
    game._path_index = 11
    game._pellets[0] = [False] * 12
    game._advance_game()
    assert game._path_index == 11
    assert not game._restore_queue
    assert not any(game._blocks)
    # While stunned, an existing recovery queue advances one cell per tick.
    game._restore_queue.extend((1, col) for col in range(7))
    game._pellets[1] = [False] * 12
    for _ in range(7):
        game._advance_game()
        assert game._path_index == 11
    assert game._caught_ticks == 0
    assert sum(game._pellets[1]) == 7
    game._ghosts.clear()
    game._advance_game()
    assert game._path[game._path_index] == (1, 0)
    assert len(game._restore_queue) == 11


def test_row_sweep_does_not_repeat_top_and_single_row_recovers():
    game = board(8, 3)
    assert [row for row, col in game._path if col == 0] == [0, 1, 2, 1]
    game = board(8, 1)
    game._ghosts.clear()
    for _ in range(8):
        game._advance_game()
    assert game._path_index == 0
    assert game._restore_queue


def test_recovery_queue_and_obstacles_stay_bounded():
    game = board()
    game._pellets[0] = [False] * 12
    game._queue_row_restore(0)
    blocks = [set(row) for row in game._blocks]
    game._queue_row_restore(0)
    assert len(game._restore_queue) == 12
    assert game._blocks == blocks
    for _ in range(100):
        game._spawn_breakable_blocks()
    assert all(len(row) <= 2 for row in game._blocks)


def test_narrow_ghost_wraps_and_giant_respawn_avoids_walls():
    game = board(2, 3)
    ghost = game._ghosts[0]
    game._ghosts = [ghost]
    ghost.row, ghost.column, ghost.direction = 0, 0, 1
    game._move_ghosts(2, 0, 2)
    assert (ghost.row, ghost.column, ghost.direction) == (1, 0, 1)
    game._blocks[0] = {0}
    game._giant_ticks = 10
    game._handle_collisions(1, 0)
    assert ghost.row == 2


def test_frightened_ghost_can_leave_a_wall_after_bouncing():
    game = board()
    ghost = game._ghosts[0]
    game._ghosts = [ghost]
    ghost.row, ghost.column, ghost.direction = 0, 3, -1
    game._blocks[0] = {2}
    game._giant_ticks = 10
    game._move_ghosts(0, 8, 12)
    assert ghost.direction == 1
    game._move_ghosts(0, 8, 12)
    assert ghost.column == 4


def test_rendering_width_and_state_survive_long_cycles_and_resizes():
    game = PacmanGame()
    for width in (1, 2, 3, 12, 55):
        for height in (1, 2, 4):
            game._ensure_board(width, height)
            assert game._caught_ticks == game._giant_ticks == 0
            for tick in range(250):
                game._step += 1
                game._advance_game()
                lines = game._render_board().split("\n")
                assert len(lines) == height
                assert all(visible_width(line) == width for line in lines)
                assert len(game._restore_queue) <= width * height
                assert all(len(row) <= 2 for row in game._blocks)
            # Force both blink phases at the final column, including 1-cell space.
            game._path_index = width - 1
            game._caught_ticks = 7
            for phase in (0, 1):
                game._step = phase
                assert all(visible_width(line) == width for line in game._render_board().split("\n"))
