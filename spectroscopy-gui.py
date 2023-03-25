""" An application for interactive treatment and analysis of spectra

This graphicaly interfaced application is conceived to use the spectroscopy
library, matplotlib and gui library to provide an easy to use application
for the treatment and analysis of spectra from different sources.
"""

import spectroscopy as spec
import spectroscopy.gui as gui

if __name__ == '__main__':
    gui.splash()
    gui.interface()
    gui.bye()
print( messages.ABOUT )
