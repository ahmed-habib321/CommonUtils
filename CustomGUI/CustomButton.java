package CustomGUI;

import java.awt.Color;
import java.awt.Cursor;
import java.awt.Graphics;

import javax.swing.JButton;

/**
 * A custom button component with a rounded rectangle shape and a configurable
 * border.
 */
public class CustomButton extends JButton {

    private int red;
    private int green;
    private int blue;
    private Color border;

    /**
     * Constructs a new CustomButton instance with default properties.
     */
    public CustomButton() {
        setProperties(Color.WHITE, new Cursor(Cursor.HAND_CURSOR));
        this.red = 70;
        this.green = 130;
        this.blue = 180;
        this.border = Color.BLUE;
    }

    /**
     * Constructs a new CustomButton instance with the specified properties.
     *
     * @param foreground the foreground color of the button
     * @param cursor the cursor type of the button
     * @param red the red component of the button's background color
     * @param green the green component of the button's background color
     * @param blue the blue component of the button's background color
     * @param border the color of the button's border
     */
    public CustomButton(Color foreground, Cursor cursor, int red, int green, int blue, Color border) {
        setProperties(foreground, cursor);
        this.red = red;
        this.green = green;
        this.blue = blue;
        this.border = border;
    }

    /**
     * Sets the properties of the button.
     *
     * @param foreground the foreground color of the button
     * @param cursor the cursor type of the button
     */
    private void setProperties(Color foreground, Cursor cursor) {
        setOpaque(false);
        setContentAreaFilled(false);
        setForeground(foreground);
        setCursor(cursor);
    }

    /**
     * Paints the component on the screen.
     *
     * @param g the Graphics object used for painting
     */
    @Override
    protected void paintComponent(Graphics g) {
        g.setColor(new Color(red, green, blue));
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
