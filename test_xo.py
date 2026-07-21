import unittest

from xo import display_board, is_draw, new_board, winner


class XOTests(unittest.TestCase):
    def test_new_board_is_empty(self) -> None:
        self.assertEqual(new_board(), [" "] * 9)

    def test_empty_positions_are_numbered(self) -> None:
        board = new_board()
        board[0] = "X"
        rendered = display_board(board)
        self.assertIn(" X | 2 | 3 ", rendered)
        self.assertIn(" 7 | 8 | 9 ", rendered)

    def test_winner_finds_rows_columns_and_diagonals(self) -> None:
        winning_boards = (
            ["X", "X", "X", " ", " ", " ", " ", " ", " "],
            ["O", " ", " ", "O", " ", " ", "O", " ", " "],
            ["X", " ", " ", " ", "X", " ", " ", " ", "X"],
        )
        expected_winners = ("X", "O", "X")

        for board, expected in zip(winning_boards, expected_winners):
            with self.subTest(board=board):
                self.assertEqual(winner(board), expected)

    def test_draw_requires_full_board_without_winner(self) -> None:
        draw = ["X", "O", "X", "X", "O", "O", "O", "X", "X"]
        self.assertTrue(is_draw(draw))
        self.assertFalse(is_draw(new_board()))


if __name__ == "__main__":
    unittest.main()
