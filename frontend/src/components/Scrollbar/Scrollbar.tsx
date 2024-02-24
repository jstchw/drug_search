import { Nav } from 'react-bootstrap';

interface ScrollbarProps {
    children: React.ReactNode;
}

const Scrollbar: React.FC<ScrollbarProps> = ({ children }) => {
    return (
        <Nav
            className="flex-column fixed-top bg-light"
        >
            {children}
        </Nav>
    );
}

export default Scrollbar;