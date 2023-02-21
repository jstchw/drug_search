import Navbar from 'react-bootstrap/Navbar';
import Container from 'react-bootstrap/Container';
import { NavLink } from 'react-router-dom';
import '../styles/main.css'
import HeaderLink from './HeaderLink'

const Header = () => {
    return (
        // Navbar containing the links to all the pages of the webapp
        <Navbar bg="dark" variant="dark" sticky="top" className="header-full-width">
            <Container>
                <Navbar.Brand as={NavLink} to="/">DrugSearch</Navbar.Brand>
                <HeaderLink to="/dashboard">Dashboard</HeaderLink>
            </Container>
        </Navbar>
    );
}

export default Header;