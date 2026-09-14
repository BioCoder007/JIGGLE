import os
import sys

from PySide6.QtGui import QIcon
from PySide6.QtWidgets import QApplication

from ui.main_window import MainWindow

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    app = QApplication(sys.argv)

    icon_path = os.path.join(HERE, "ui", "assets", "icon.png")
    if os.path.exists(icon_path):
        app.setWindowIcon(QIcon(icon_path))

    win = MainWindow()
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
