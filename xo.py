"""A simple two-player command-line XO (tic-tac-toe) game."""

from __future__ import annotations

WINNING_LINES = (
    (0, 1, 2),
    (3, 4, 5),
    (6, 7, 8),
    (0, 3, 6),
    (1, 4, 7),
    (2, 5, 8),
    (0, 4, 8),
    (2, 4, 6),
)


def new_board() -> list[str]:
    """Return an empty board."""
    return [" "] * 9


def display_board(board: list[str]) -> str:
    """Format the board, showing position numbers for empty squares."""
    cells = [value if value != " " else str(index + 1) for index, value in enumerate(board)]
    return (
        f" {cells[0]} | {cells[1]} | {cells[2]} \n"
        "---+---+---\n"
        f" {cells[3]} | {cells[4]} | {cells[5]} \n"
        "---+---+---\n"
        f" {cells[6]} | {cells[7]} | {cells[8]} "
    )


def winner(board: list[str]) -> str | None:
    """Return the winning mark, or None if nobody has won."""
    for first, second, third in WINNING_LINES:
        if board[first] != " " and board[first] == board[second] == board[third]:
            return board[first]
    return None


def is_draw(board: list[str]) -> bool:
    """Return whether the board is full without a winner."""
    return winner(board) is None and all(cell != " " for cell in board)


def read_move(board: list[str], player: str) -> int:
    """Prompt until the player chooses an available square."""
    while True:
        try:
            choice = input(f"Player {player}, choose a square (1-9): ").strip()
        except (EOFError, KeyboardInterrupt):
            print("\nGame ended.")
            raise SystemExit(0)

        if not choice.isdigit() or not 1 <= int(choice) <= 9:
            print("Please enter a number from 1 to 9.")
            continue

        position = int(choice) - 1
        if board[position] != " ":
            print("That square is already taken.")
            continue

        return position


def play_game() -> None:
    """Run one game of XO."""
    board = new_board()
    player = "X"

    print("XO (Tic-Tac-Toe)")
    print("Choose a numbered square to place your mark.\n")

    while True:
        print(display_board(board))
        position = read_move(board, player)
        board[position] = player

        winning_player = winner(board)
        if winning_player:
            print(f"\n{display_board(board)}")
            print(f"\nPlayer {winning_player} wins!")
            return

        if is_draw(board):
            print(f"\n{display_board(board)}")
            print("\nIt's a draw!")
            return

        player = "O" if player == "X" else "X"
        print()


if __name__ == "__main__":
    play_game()
