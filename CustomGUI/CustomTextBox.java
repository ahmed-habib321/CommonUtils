package CustomGUI;

import java.awt.Color;
import java.awt.Graphics;
import javax.swing.JTextField;

public class CustomTextBox extends JTextField {

    private Color border;

    /**
     * Constructs a new CustomPasswordBox instance with default background color
     * and red border.
     */
    public CustomTextBox() {
        setProperties();
        this.border = Color.RED;
    }

    /**
     * Constructs a new CustomTextBox instance with the specified RGB values and
     * border color.
     *
     * @param border the color of the border
     */
    public CustomTextBox(Color border) {
        setProperties();
        this.border = border;
    }

    /**
     * Sets the properties of the Text field.
     */
    private void setProperties() {
        setOpaque(false);
    }

    /**
     * Paints the component on the screen.
     *
     * @param g the Graphics object used for painting
     */
    @Override
    protected void paintComponent(Graphics g) {
        g.setColor(getBackground());
        g.fillRoundRect(0, 0, getWidth() - 1, getHeight() - 1, getHeight() - 1, getHeight() - 1);
        super.paintComponent(g);
    }

    /**
     * Paints the border of the component.
     *
     * @param g the Graphics object used for painting
     */
    @Override
    protected void paintBorder(Graphics g) {
        g.setColor(border);
        g.drawRoundRect(0, 0, getWidth() - 1, getHeight() - 1, getHeight() - 1, getHeight() - 1);
    }
}
