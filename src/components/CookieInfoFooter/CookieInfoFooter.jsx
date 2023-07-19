import {Button, Col, Row} from "react-bootstrap";
import {useEffect, useState} from "react";
import './CookieInfoFooter.css'
import {CSSTransition} from "react-transition-group";

export const CookieInfoFooter = () => {
    const [show, setShow] = useState(false)

    useEffect(() => {
        const dismissed = localStorage.getItem('cookie-info-dismissed')
        if (!dismissed) {
            setShow(true)
        }
    }, [])

    const handleDismissed = () => {
        localStorage.setItem('cookie-info-dismissed', true)
        setShow(false)
    }

    return (
        <CSSTransition
            in={show}
            timeout={300}
            classNames="fade"
            unmountOnExit
        >
            { show ? (
            <div className={'cookie-info-footer fixed-bottom d-flex justify-content-center m-4 text-center'}>
                <Row>
                    <Col md="auto" className={'p-3'}>
                        We use cookies for functionality, no personal data is collected.
                        <Button className={'mx-2'} onClick={handleDismissed}>Got it!</Button>
                    </Col>
                </Row>
            </div> ) : <></> }
        </CSSTransition>
    )
}